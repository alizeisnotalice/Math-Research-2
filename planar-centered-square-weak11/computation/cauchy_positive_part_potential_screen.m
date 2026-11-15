function cauchy_positive_part_potential_screen(output_root)
% Numeric screening only. Full positive-part Cauchy potential, including infinite axial tails.
assert(startsWith(output_root,'/Volumes/data/'));if ~isfolder(output_root),mkdir(output_root);end
weights=[0,.25,.5,.75,1];Rvals=1.05:.01:4;tvals=(0:64)/128;
rows=struct([]);started=tic;
for k=1:numel(weights)
    theta=weights(k);
    integrand=@(t,y) 8*max(potential(y+t./(1-t),y,theta),0)./(1-t).^2;
    I=integral2(integrand,0,1,0,1,'Method','tiled','AbsTol',1e-7,'RelTol',1e-6);
    fprintf('integral completed theta=%.4g mass=%.15g\n',theta,I);
    scans=struct([]);cursor=0;
    for R=Rvals
        vals=[];
        for radial=[.125,.25,.5,.75,1]
            vals=[vals,potential(radial*(1-tvals)/R,radial*tvals/R,theta)]; %#ok<AGROW>
        end
        % Include exact diagonal-coordinate crossing events when they occur inside D.
        for t=tvals
            for radial=[R/(2*(1-t)),R/(2*max(t,realmin))]
                if radial>0&&radial<1,vals(end+1)=potential(radial*(1-t)/R,radial*t/R,theta);end %#ok<AGROW>
            end
        end
        gm=min(vals);if gm<=0||~isfinite(gm),continue;end
        cursor=cursor+1;scans(cursor).R=R;scans(cursor).minimum_sample_Wunit=gm;scans(cursor).sample_cost=R^2*I/(2*gm);
    end
    assert(~isempty(scans));[~,j]=min([scans.sample_cost]);
    rows(k).axis_weight=theta;rows(k).unit_positive_mass=I;rows(k).best_sample=scans(j);rows(k).R_scan=scans;
    fprintf('theta=%.4g mass=%.15g bestR=%.5g sampledCost=%.15g\n',theta,I,scans(j).R,scans(j).sample_cost);
end
out=struct('status','completed_numerical_positive_part_potential_screen','matlab_version',version,'weights',weights,'model','U=1/l1; nu=theta*four_axis_atoms+(1-theta)*four_diagonal_atoms at common l1 radius R; K=(U-U*nu)_+','integral','8*int_y=0^1 int_t=0^1 Kunit(y+t/(1-t),y)/(1-t)^2 dt dy; infinite tails retained','integral_AbsTol',1e-7,'integral_RelTol',1e-6,'parameter_samples','theta 0,.25,.5,.75,1; R1.05:.01:4; 65 boundary directions times five radial layers and diagonal coordinate crossings','models',rows,'elapsed_seconds',toc(started),'scope','Float quadrature and finite D samples only. No integral enclosure or continuous D majorization certificate; no new general C2 upper.');
fid=fopen(fullfile(output_root,'positive-part-potential-screen.json'),'w');assert(fid>=0);fwrite(fid,jsonencode(out,PrettyPrint=true));fwrite(fid,newline);fclose(fid);
end
function w=potential(x,y,theta)
if theta==1,w=(pair(x,y,1,0)+pair(x,y,0,1))/2;return;end
if theta==0,w=(pair(x,y,.5,.5)+pair(x,y,.5,-.5))/2;return;end
wa=(pair(x,y,1,0)+pair(x,y,0,1))/2;
wd=(pair(x,y,.5,.5)+pair(x,y,.5,-.5))/2;
w=theta*wa+(1-theta)*wd;
end
function w=pair(x,y,a,b)
r=x+y;dm=delta(x,a)+delta(y,b);dp=delta(x,-a)+delta(y,-b);
w=.5*(dm./(r.*(r+dm))+dp./(r.*(r+dp)));
end
function d=delta(u,c)
if c>=0,d=c-2*min(u,c);else,d=(-c)+zeros(size(u));end
end
