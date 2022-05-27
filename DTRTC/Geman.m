function [gamma1,gamma2]=Geman(t1,t2,geman)

gamma1=geman/(t1+geman)^2; gamma2=geman/(t2+geman)^2;

end