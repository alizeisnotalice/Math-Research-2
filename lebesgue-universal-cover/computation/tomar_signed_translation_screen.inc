        // The signed area of this fixed material chain is separately affine
        // in all eight translations when all three angles are fixed.
        // Prove monotonicity first, then only enumerate undecided endpoints.
        if (twice-error < goal) {
            L tx[5]{},ty[5]{};
            int type[8]{},body[8]{},free_count=0;
            for (int k=1;k<=4;++k) {
                W variation_x=0,variation_y=0;
                for (int j=1;j<=4;++j) {
                    int pair=std::abs(translation_pairs[std::min(k,j)][std::max(k,j)]);
                    variation_x+=W(pair)*hy[j];
                    variation_y+=W(pair)*hx[j];
                }
                for (int j=2;j<=4;++j) {
                    const V v=mixed[k][j];
                    variation_x+=up(W(std::abs(v.x))*theta[j]+W(std::abs(v.y))*cos_width[j],Q);
                    variation_y+=up(W(std::abs(v.y))*theta[j]+W(std::abs(v.x))*cos_width[j],Q);
                }
                auto fix=[&](int which,L gradient,W variation,L width,L& out) {
                    if (gradient>=variation) out=-width;
                    else if (gradient<=-variation) out=width;
                    else {type[free_count]=which;body[free_count++]=k;}
                };
                fix(0,g[k].y,variation_x,hx[k],tx[k]);
                fix(1,-g[k].x,variation_y,hy[k],ty[k]);
            }
            if (free_count <= 4) {
                W worst=0;
                for (int bits=0;bits<(1<<free_count);++bits) {
                    for(int i=0;i<free_count;++i) {
                        const int k=body[i],sign=((bits>>i)&1)?1:-1;
                        if(type[i])ty[k]=sign*hy[k];else tx[k]=sign*hx[k];
                    }
                    W value=twice-rotation_remainder,sg[5]{},cs[5]{};
                    for(int j=2;j<=4;++j) {sg[j]=grot[j];cs[j]=cosine[j];}
                    for(int k=1;k<=4;++k) {
                        value+=W(tx[k])*g[k].y-W(ty[k])*g[k].x;
                        for(int j=k+1;j<=4;++j)
                            value+=W(translation_pairs[k][j])*det(tx[k],ty[k],tx[j],ty[j]);
                        for(int j=2;j<=4;++j) {
                            auto v=mixed[k][j];
                            sg[j]+=W(tx[k])*v.x+W(ty[k])*v.y;
                            cs[j]+=det(tx[k],ty[k],v.x,v.y);
                        }
                    }
                    for(int j=2;j<=4;++j) {
                        value-=up(ab(sg[j])*(sg[j]>=0?negative_angle[j]:positive_angle[j]),Q);
                        if(cs[j]>0)value-=up(cs[j]*cos_width[j],Q);
                    }
                    if(bits==0 || value<worst)worst=value;
                }
                if(worst>twice-error)error=twice-worst;
            }
        }
