function [o_noise_sigma, o_offset] = f_estimate_noise_and_offset(i_image)

%A smartass way to find out the hardware background offset and noise sigma
%
%Please use with extreme care - this function makes A LOT of assumpitons
%about input (3D, mode is background offset, Gaussian distribution in the
%bg, maybe more)
%{
 A crazy-lazy hack to find the SD (sigma) of the noise- exchanging programmer time for machine time:
    Assuming the noise is gaussian with the mean of t_hwbg,
    let's take the left part of the this distribution, found with
    A(A<t_hwbg), and reflect it to the right, and then use the standard SD
function on the resulting stuff
%}
o_offset=mode(mode(mode(i_image))); % this might be rewritten as a while numel()>1 loop for higher dims
t_left=i_image(i_image<o_offset);
%hist([t_hwbg*2-left;t_left]);figure(gcf); %for test
o_noise_sigma=std([o_offset*2-t_left;t_left]);