function v=DSS_Aux_LogUniform(x,a,b)
v = (x>a & x<b).*1./(x.*(log10(b) - log10(a)));