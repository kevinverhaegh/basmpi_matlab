function Y=DSS_Aux_LogUniform_Sample(N,a,b)
Iter = 1e5;
if numel(N)==1
    N = rand(N,1);
end
x=logspace(log10(a),log10(b),Iter);
v=DSS_Aux_LogUniform(x,a,b);
vc = cumtrapz(x,v);
vc = vc./max(vc);
Y = interp1(vc,x,N);