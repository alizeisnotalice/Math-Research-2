function cauchy_mixedtime_source_lp(output_root)
% R60 exploratory joint source/time LP. Numerical cell coefficients, not a certificate.
assert(startsWith(output_root,'/Volumes/data/'));if ~isfolder(output_root),mkdir(output_root);end
started=tic;radii=[1.25,1.5,2,3,4,8,16];shapes=[0,.25,.5];times=[0,.1,.4,1];orbits=struct('R',{},'shape',{},'time',{},'points',{});J=0;
for R=radii,for shape=shapes
 pts=[];a=R*(1-shape);b=R*shape;for sx=[-1,1],for sy=[-1,1],pts=[pts;sx*a,sy*b;sx*b,sy*a];end,end %#ok<AGROW>
 pts=unique(pts,'rows');for tt=times,J=J+1;orbits(J)=struct('R',R,'shape',shape,'time',tt,'points',pts);end
end,end
for tt=[.1,.4,1],J=J+1;orbits(J)=struct('R',0,'shape',0,'time',tt,'points',[0,0]);end
L=max(radii);coords=[];for j=1:J,coords=[coords;abs(orbits(j).points(:))];end %#ok<AGROW>
yn=unique([0:1/8:4,4.5:.5:8,9:L,coords']);xn=L:2:4*L;ny=numel(yn)-1;na=numel(xn)-1;ncentral=ny^2;narm=(na+1)*ny;ncell=ncentral+narm;
[X,Y]=ndgrid(yn,yn);[XA,YA]=ndgrid(xn,yn);mult=[2*ones(ncentral,1);4*ones(narm,1)];C=zeros(ncell,J);[gn,gw]=gauss(24);
for j=1:J
 pts=orbits(j).points;tt=orbits(j).time;n=size(pts,1);P=primitive(X,Y);PA=primitive(XA,YA);M=-sxlogx(yn);
 for k=1:n,P=P-primitive(X-pts(k,1),Y-pts(k,2))/n;PA=PA-primitive(XA-pts(k,1),YA-pts(k,2))/n;M=M+sxlogx(yn-pts(k,2))/n;end
 cc=diff(diff(P,1,1),1,2);ca=diff(diff(PA,1,1),1,2);ct=diff(M-PA(end,:));
 if tt>0
  for k=1:n
   z=pts(k,:);
   for q=1:numel(gn)
    s=tt*gn(q);fac=2*tt*gw(q)/(pi*n);dx=angle_diff(yn-z(1),s);dy=angle_diff(yn-z(2),s);da=angle_diff(xn-z(1),s);
    cc=cc+fac*(dx(:)*dy(:)');ca=ca+fac*(da(:)*dy(:)');ct=ct+fac*atan2(s,xn(end)-z(1))*dy;
   end
  end
 end
 C(:,j)=[cc(:);ca(:);ct(:)];
end
% The remaining full-plane corner has exact signed half mass pi*t - covered half mass.
comp=pi*[orbits.time]-mult'*C;
[dx,dy]=diamond_points(32,128);WD=zeros(numel(dx),J);for j=1:J,WD(:,j)=basis(dx,dy,orbits(j));end
Amat=[sparse(-WD),sparse(numel(dx),ncell+1);sparse(C),-speye(ncell),sparse(ncell,1);sparse(comp),sparse(1,ncell),-1];
bvec=[-ones(numel(dx),1);zeros(ncell+1,1)];f=[zeros(J,1);mult;1];opts=optimoptions('linprog','Display','iter','ConstraintTolerance',1e-8,'OptimalityTolerance',1e-8);
[a,fval,ef,lpout]=linprog(f,Amat,bvec,[],[],zeros(J+ncell+1,1),[],opts);assert(ef>0);coef=max(a(1:J),0);
[fx,fy]=diamond_points(128,512);fv=total(fx,fy,coef,orbits);[dm,ix]=min(fv);arg=[fx(ix),fy(ix)];
o=optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',200,'OptimalityTolerance',1e-10,'StepTolerance',1e-12);[~,ix]=sort(fv);local=[];
for q=1:5
 [zz,mz,flag]=fmincon(@(z)total(z(1),z(2),coef,orbits),[fx(ix(q)),fy(ix(q))],[-1,1;1,1;-1,-1],[0;1;-1e-6],[],[],[0,0],[1,1],[],o);local=[local;zz,mz,flag];if mz<dm,dm=mz;arg=zz;end %#ok<AGROW>
end
assert(dm>0);coef=coef/dm;
fcore=@(u,y)4*L*max(total(y+L*u./(1-u),y,coef,orbits),0)./(1-u).^2;
core=integral2(fcore,0,1,0,L,'Method','tiled','AbsTol',2e-5,'RelTol',2e-4);
ftail=@(u,v)4*L^2*max(total(L+L*v./(1-v)+L*u./(1-u),L+L*v./(1-v),coef,orbits),0)./((1-u).^2.*(1-v).^2);
tail=integral2(ftail,0,1,0,1,'Method','tiled','AbsTol',2e-5,'RelTol',2e-4);
for j=1:J,orbits(j).coefficient=coef(j);end
out=struct('status','completed_numerical_joint_source_time_LP','matlab_version',version,'threads',1,'radii',radii,'shapes',shapes,'times',times,'source_count',J,'cell_count',ncell,'quadrature_order_for_Qcells',24,'LP_exitflag',ef,'LP_output',lpout,'LP_objective_float_cell_approximation',fval,'D_coarse_count',numel(dx),'D_fine_count',numel(fx),'minimum_before_rescale',dm,'minimum_point',arg,'local_optimizers',local,'actual_numeric_half_integral',core+tail,'core_half_integral',core,'corner_half_integral',tail,'total_mass',sum(coef),'weighted_time',coef'*[orbits.time]','signed_half_mass',pi*coef'*[orbits.time]','orbits',orbits,'y_nodes',yn,'arm_x_nodes',xn,'elapsed_seconds',toc(started),'scope','Cell objective is a formal Jensen lower model; its coefficients use undirected Gauss24. It is not certified even as lower bound. Full-D finite samples/local search and whole-plane quadrature are exploratory only; no general C2 upper.');
fid=fopen(fullfile(output_root,'mixedtime-source-lp.json'),'w');assert(fid>=0);fwrite(fid,jsonencode(out,PrettyPrint=true));fwrite(fid,newline);fclose(fid);
fprintf('R60 J=%d cells=%d LP=%.15g fullnumeric=%.15g fineMin=%.15g active=%d elapsed=%.4g\n',J,ncell,fval,core+tail,dm,nnz(coef),toc(started));
end
function v=basis(x,y,o)
n=size(o.points,1);v=potential(x,y,o.time,o.points,ones(n,1)/n);
end
function v=total(x,y,a,orbits)
v=zeros(size(x));for j=1:numel(a),if a(j)>0,v=v+a(j)*basis(x,y,orbits(j));end,end
end
function d=angle_diff(edges,s)
d=atan2(s*diff(edges),s*s+edges(1:end-1).*edges(2:end));
end
function [x,w]=gauss(n)
b=(1:n-1)./sqrt(4*(1:n-1).^2-1);[V,D]=eig(diag(b,1)+diag(b,-1));[x,ix]=sort(diag(D));x=(x+1)/2;w=V(1,ix).^2;
end
function h=primitive(x,y)
a=abs(x);b=abs(y);h=sign(x).*sign(y).*(xlogx(a+b)-xlogx(a)-xlogx(b));
end
function f=xlogx(x)
f=zeros(size(x));i=x>0;f(i)=x(i).*log(x(i));
end
function f=sxlogx(x)
f=zeros(size(x));i=x~=0;f(i)=x(i).*log(abs(x(i)));
end
function [x,y]=diamond_points(N,Nb)
x=[];y=[];for i=0:N,for j=0:i,if i+j<=N&&i+j>0,x(end+1,1)=i/N;y(end+1,1)=j/N;end,end,end %#ok<AGROW>
h=(0:Nb)'/(2*Nb);p=unique([x,y;1-h,h],'rows');x=p(:,1);y=p(:,2);
end
function w=potential(x,y,t,pts,mass)
r=x+y;w=zeros(size(r));
if t==0
 first=zeros(size(r));for j=1:numel(mass)
 d=delta(x,pts(j,1))+delta(y,pts(j,2));
 if pts(j,1)>0,first=first+2*mass(j)*max(pts(j,1)-x,0);end
 if pts(j,2)>0,first=first+2*mass(j)*max(pts(j,2)-y,0);end
 w=w-mass(j)*d.^2./(r.^2.*(r+d));end
 w=w+first./r.^2;return;
end
near=r<64;far=~near;
if any(near(:))
 xx=x(near);yy=y(near);v=sum(mass)./(xx+yy);
 for j=1:numel(mass),v=v-mass(j)*smoothU(abs(xx-pts(j,1)),abs(yy-pts(j,2)),t);end
 w(near)=v;
end
if any(far(:))
 xx=x(far);yy=y(far);rr=xx+yy;v=zeros(size(rr));first=zeros(size(rr));
 for j=1:numel(mass)
  a=abs(xx-pts(j,1));b=abs(yy-pts(j,2));d=delta(xx,pts(j,1))+delta(yy,pts(j,2));
  if pts(j,1)>0,first=first+2*mass(j)*max(pts(j,1)-xx,0);end
  if pts(j,2)>0,first=first+2*mass(j)*max(pts(j,2)-yy,0);end
  v=v+mass(j)*(-d.^2./(rr.^2.*(rr+d))+shortQ(a,b,t));
 end
 w(far)=v+first./rr.^2;
end
end
function v=smoothU(a,b,t)
lo=min(a,b);hi=max(a,b);h=hi-lo;den=t*t+lo.*hi;z=t*h./den;ratio=ones(size(z));ii=z~=0;ratio(ii)=atan(z(ii))./z(ii);
v=(2/pi)*(atan(hi/t)+lo*t./den.*ratio)./(lo+hi);zz=lo+hi==0;v(zz)=2/(pi*t);
end
function q=shortQ(a,b,t)
lo=min(a,b);hi=max(a,b);q=zeros(size(lo));ii=lo<=8*t;
if any(ii(:))
 aa=lo(ii);bb=hi(ii);h=bb-aa;den=t*t+aa.*bb;z=t*h./den;ratio=ones(size(z));jj=z~=0;ratio(jj)=atan(z(jj))./z(jj);
 q(ii)=(2/pi)*(atan(t./bb)-aa*t./den.*ratio)./(aa+bb);
end
jj=~ii;
if any(jj(:))
 a2=lo(jj).^2;b2=hi(jj).^2;
 nodes=[.0198550717512319,.1016667612931866,.2372337950418355,.4082826787521751,.5917173212478249,.7627662049581645,.8983332387068134,.9801449282487681];
 wt=[.0506142681451881,.1111905172266872,.1568533229389437,.181341891689181,.181341891689181,.1568533229389437,.1111905172266872,.0506142681451881];
 v=zeros(size(a2));for n=1:8,s2=(t*nodes(n))^2;v=v+wt(n)*s2./((s2+a2).*(s2+b2));end;q(jj)=(2*t/pi)*v;
end
end
function d=delta(u,c)
if c>=0,d=c-2*min(u,c);else,d=(-c)+zeros(size(u));end
end
