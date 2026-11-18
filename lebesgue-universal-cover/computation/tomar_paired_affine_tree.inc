    // Affine coefficients use relative box coordinates in [0,1]^11.
    using Plane=std::array<L,12>;
    struct DirectionPlanes {std::array<Plane,5> planes{};int size=0;L midpoint=0;};
    std::vector<DirectionPlanes> direction_planes;
    std::vector<L> direction_stamp;
    Plane affine_area{};
    L stamp=0,paired_boxes=0,paired_cells=0,paired_tries=0,paired_rows=0,paired_tangents=0;
    bool enable_pairing=false;

    static L lower_round(W a,L b){return -up_signed(-a,b);}
    static W min_plane(const Plane& f){
        W out=f[0];
        for(int j=1;j<12;++j)out+=std::min<L>(0,f[j]);
        return out;
    }
    static W midpoint_twice(const Plane& f){
        W out=W(2)*f[0];
        for(int j=1;j<12;++j)out+=f[j];
        return out;
    }
    void prepare_pairing(){
        if(direction_planes.empty()){
            direction_planes.resize(n);
            direction_stamp.resize(n);
        }
        constexpr L den=1024*Fine;
        affine_area[0]=t.affine_center;
        for(int j=0;j<11;++j){
            L a=lower_round(W(t.affine_coeff[j])*(1024*pair_lo[j]-t.affine_anchor[j]),den);
            L b=lower_round(W(t.affine_coeff[j])*(1024*pair_hi[j]-t.affine_anchor[j]),den);
            affine_area[0]+=a;
            affine_area[j+1]=b-a;
        }
        ++stamp;
    }
    Box pair_lo{},pair_hi{};
    const DirectionPlanes& body_planes(int j){
        auto& result=direction_planes[j];
        if(direction_stamp[j]==stamp)return result;
        direction_stamp[j]=stamp;
        const L x=directions[2*j],y=directions[2*j+1];
        constexpr L den=(Fine/Q)*U;
        std::array<Plane,5> raw{};
        raw[0][0]=disk[j];
        for(int i=0;i<4;++i){
            auto& f=raw[i+1];
            int k=3*i;
            const L x0=up_signed(W(x)*pair_lo[k],den), x1=up_signed(W(x)*pair_hi[k],den);
            const L y0=up_signed(W(y)*pair_lo[k+1],den), y1=up_signed(W(y)*pair_hi[k+1],den);
            f[0]=x0+y0;
            f[k+1]=x1-x0;
            f[k+2]=y1-y0;
            const L h=endpoints[2*i*n+j],p=endpoints[(2*i+1)*n+j];
            if(i){
                --k;
                auto rotate=[&](L a){
                    return up_signed(W(h)*(root_hi[k]-a)+W(p)*(a-root_lo[k]),root_hi[k]-root_lo[k]);
                };
                const L h0=rotate(pair_lo[k]),h1=rotate(pair_hi[k]);
                f[0]+=h0;
                f[k+1]=h1-h0;
            }else f[0]+=h;
        }
        result.size=0;
        result.midpoint=std::numeric_limits<L>::min();
        for(unsigned i=0;i<raw.size();++i){
            bool keep=true;
            for(unsigned k=0;k<raw.size();++k){
                if(k==i)continue;
                if(raw[k]==raw[i]){
                    if(k<i){keep=false;break;}
                    continue;
                }
                W margin=W(raw[k][0])-raw[i][0];
                for(int a=1;a<12;++a)margin+=std::min<W>(0,W(raw[k][a])-raw[i][a]);
                if(margin>=0){keep=false;break;}
            }
            if(keep){
                result.planes[result.size++]=raw[i];
                result.midpoint=std::max(result.midpoint,L(midpoint_twice(raw[i])));
            }
        }
        if(!result.size)std::abort();
        return result;
    }
    bool tangent_positive(const DirectionPlanes& a,const DirectionPlanes& b,
                          const DirectionPlanes& c,const L* q,L slope,L penalty) const {
        std::array<W,12> base{};
        const L w0=q[6];
        base[0]=W(w0)*(W(affine_area[0])-penalty-t.goal+W(slope)*q[7]);
        for(int j=1;j<12;++j)base[j]=W(w0)*affine_area[j];
        std::array<W,3> weight={W(slope)*q[3],W(slope)*q[4],W(slope)*q[5]};
        for(int i=0;i<a.size;++i){
            const auto& f=a.planes[i];
            for(int j=0;j<b.size;++j){
                const auto& g=b.planes[j];
                std::array<W,12> first{};
                for(int d=0;d<12;++d)first[d]=base[d]-weight[0]*f[d]-weight[1]*g[d];
                for(int k=0;k<c.size;++k){
                    const auto& h=c.planes[k];
                    W lower=first[0]-weight[2]*h[0];
                    for(int d=1;d<12;++d)lower+=std::min<W>(0,first[d]-weight[2]*h[d]);
                    if(lower<=0)return false;
                }
            }
        }
        return true;
    }
    template<class Support,class Threshold>
    bool paired_cell(int i,L e,Support support_upper,Threshold threshold){
        ++paired_cells;
        const W area2=midpoint_twice(affine_area);
        for(L k=pointers[i];k<pointers[i+1];++k){
            ++paired_rows;
            const L* q=triples+8*k;
            const auto& a=body_planes(q[0]);
            const auto& b=body_planes(q[1]);
            const auto& c=body_planes(q[2]);
            const W m2=W(2)*q[6]*q[7]-W(q[3])*a.midpoint-W(q[4])*b.midpoint-W(q[5])*c.midpoint;
            if(m2<=0)continue;
            const L m=L(m2/(W(2)*q[6]));
            if(!m)continue;
            W norm2=0;
            L stops[3];
            for(int j=0;j<3;++j){
                const L d=threshold(q[j])-support_upper(q[j]);
                stops[j]=up_signed(W(q[j+3])*d,q[6]);
                norm2+=W(stops[j])*stops[j];
            }
            auto pair=overlap[k];
            if(pair[0]>=0)norm2+=W(2)*stops[pair[0]]*stops[pair[1]];
            if(norm2<=0)std::abort();
            if(area2*norm2+W(2)*e*(m+1)*(m+1)<=W(2)*t.goal*norm2)continue;
            for(int factor:{4,3,5,2,6}){
                W candidate=W(factor)*e*m/(2*norm2);
                if(candidate>std::numeric_limits<L>::max())continue;
                L slope=L(candidate);
                if(!slope)continue;
                L penalty=up_signed(W(slope)*slope*norm2,4*e);
                ++paired_tangents;
                if(tangent_positive(a,b,c,q,slope,penalty))return true;
            }
        }
        return false;
    }
