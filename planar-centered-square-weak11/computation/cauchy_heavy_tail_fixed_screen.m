function cauchy_heavy_tail_fixed_screen(input_path,output_root)
% Floating decision screen. Full unbounded wedge, no tail silently dropped.
assert(startsWith(output_root,'/Volumes/data/'));if ~isfolder(output_root),mkdir(output_root);end
started=tic;inp=jsondecode(fileread(input_path));orbits=inp.orbits;pts=[];mass=[];
for j=1:numel(orbits),z=orbits(j).points;pts=[pts;z];mass=[mass;repmat(orbits(j).coefficient/size(z,1),size(z,1),1)];end %#ok<AGROW>
assert(max(sum(abs(pts),2))<64,'Far stable branch requires source rho<64');
L=max(abs(pts(:)));times=[.025,.1,.3,1];rows=struct([]);
[x,y]=diamond_points(128,512);opts=optimoptions('fmincon','Display','off','Algorithm','sqp','OptimalityTolerance',1e-10,'StepTolerance',1e-12,'MaxIterations',150);
for k=1:numel(times)
 t=times(k);v=potential(x,y,t,pts,mass);[m,ii]=min(v);loc=[x(ii),y(ii)];
 [~,ix]=sort(v);optim=[];
 for j=1:5
  st=[x(ix(j)),y(ix(j))];[zz,mz,ef]=fmincon(@(z)potential(z(1),z(2),t,pts,mass),st,[-1,1;1,1;-1,-1],[0;1;-1e-6],[],[],[0,0],[1,1],[],opts);
  optim=[optim;zz,mz,ef];if mz<m,m=mz;loc=zz;end %#ok<AGROW>
 end
 % A full quadrant is eight identical x>=y>=0 wedges. Half integral factor is 4.
 fcore=@(u,yy)4*L*max(potential(yy+L*u./(1-u),yy,t,pts,mass),0)./(1-u).^2;
 core=integral2(fcore,0,1,0,L,'Method','tiled','AbsTol',2e-5,'RelTol',2e-4);
 ftail=@(u,v)4*L^2*max(potential(L+L*v./(1-v)+L*u./(1-u),L+L*v./(1-v),t,pts,mass),0)./((1-u).^2.*(1-v).^2);
 tail=integral2(ftail,0,1,0,1,'Method','tiled','AbsTol',2e-5,'RelTol',2e-4);
 rows(k).time=t;rows(k).D_sample_count=numel(x);rows(k).minimum_sample_or_local=m;rows(k).minimum_point=loc;rows(k).local_optimizers=optim;rows(k).core_half_integral=core;rows(k).unbounded_corner_half_integral=tail;rows(k).numerical_cost=(core+tail)/m;rows(k).signed_integral_exact_formula_value=2*pi*t*sum(mass);
 fprintf('t=%.8g m=%.15g core=%.15g corner=%.15g cost=%.15g\n',t,m,core,tail,rows(k).numerical_cost);
end
out=struct('status','completed_numeric_fixed_weight_heavy_tail_screen','matlab_version',version,'threads',1,'input_path',input_path,'total_source_mass',sum(mass),'point_count',size(pts,1),'times',times,'models',rows,'quadrature_AbsTol',2e-5,'quadrature_RelTol',2e-4,'integral_domain','4 times x>=y>=0; split at y=L; x=y+Lu/(1-u), infinite y=L+Lv/(1-v)','stable_evaluation','Direct V near finite sources; otherwise compact W via noncancelling reflection-paired first moment minus exact quadratic remainder, plus positive Qt. Qt uses exact formula when min(a,b)<=8t and 8-point Gauss on [0,t] otherwise. No interval guarantee.','elapsed_seconds',toc(started),'scope','Floating fixed weight screening; finite full-dimensional D grid plus local optimization does not certify essential infimum. Full unbounded transformed quadrature is not a rigorous integral enclosure. No new general upper bound.');
fid=fopen(fullfile(output_root,'screen.json'),'w');assert(fid>=0);fwrite(fid,jsonencode(out,PrettyPrint=true));fwrite(fid,newline);fclose(fid);
end
function [x,y]=diamond_points(N,Nb)
x=[];y=[];for i=0:N,for j=0:i,if i+j<=N&&i+j>0,x(end+1,1)=i/N;y(end+1,1)=j/N;end,end,end %#ok<AGROW>
h=(0:Nb)'/(2*Nb);p=unique([x,y;1-h,h],'rows');x=p(:,1);y=p(:,2);
end
function w=potential(x,y,t,pts,mass)
r=x+y;w=zeros(size(r));near=r<64;far=~near;
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
