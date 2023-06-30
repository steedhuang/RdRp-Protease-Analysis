function [outputArg1] = Langlands(Wuhan,z)
% Automorphic function
Seq = (Wuhan(:,3)+z*Wuhan(:,4))./(Wuhan(:,1)+z*Wuhan(:,2));
% Average of mass enthalpy to charge potential 
avg = mean(Seq);
% Fractal moment for 11th string world
outputArg1 = sum((Seq - avg).^(1/11))/(length(Seq)-1);
end

