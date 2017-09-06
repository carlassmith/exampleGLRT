function [pMapP, pMapFA,L] = doBgIRatio(Npixels,sigma,I0,bgc,Nsims)   
    iterations = 30;
    
    dataH0=noise(ones(Npixels,Npixels,Nsims)*(bgc),'poisson');
    [~, ~, L]=gpuGaussMLEv3(permute(single(dataH0),[1 2 3]),single(sigma),iterations,1,false);
    PFA0 = L(3,:);
    PD0 = L(4,:);    

    coords=repmat((Npixels-1)/2,[Nsims,2]);
    coords=coords-0.5+1.*rand(Nsims,2);
    
    M1 = finiteGaussPSFerf(Npixels,sigma,I0,bgc,coords);
    dataH1=noise(M1,'poisson');
    [~, ~, L]=gpuGaussMLEv3(permute(single(dataH1),[1 2 3]),single(sigma),iterations,1,false);
    PFA1 = L(3,:);
    PD1 = L(4,:);    
    
    pMapP = [PD0 PD1]'; 
    pMapFA = [PFA0 PFA1]';
end
