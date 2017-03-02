function [out]=f_set_negatives_to_zero(in)

if any(in<0);
    fprintf(2,'negatives set to 0\n');
    out=in;
    out(in<0)=0;
else
    out=in;
end
end