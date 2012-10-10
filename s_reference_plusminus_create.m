
%% unmix the bleaching spectrum using all excitation channels that don't excite the donor and all
%non-fret spectra
    % load non-fret spectra
        % load the donor spectrum
            S_donor=load('Ref_CFP','S');
            S_donor=S_donor.S;
        % load the acceptor spectrum
            S_acceptor=load('Ref_YFP');
            S_acceptor=S_acceptor.Sa;
        % load the autofluorescence spectrum
            S_auto=load('Ref_auto');
            S_auto=S_auto.S;
        % check that indeces match
            if S_donor.wl~=S_acceptor.wl; error('Acceptor and donor indices mismatch');end           
            if S_donor.emch~=S_acceptor.emch;error('Acceptor and donor indices mismatch');end  
            if S_donor.wl~=S_auto.wl;error('Autofluorescence and donor indices mismatch');end
            if S_donor.emch~=S_auto.emch;error('Autofluorescence and donor indices mismatch');end
            %t_concentrations=nonnegative_unmix(,)
    % load the bleach spectrum 
         load('temp_spectra','temp_bleach_spectrum');
         S_bleach=temp_bleach_spectrum';
    %form the S matrix for nonnegative_unmix
        X=S_bleach(:,3:6);
        %%
    %form the X matrix for nonnegative_unmix
            S=[S_donor.spectrum(3:6); S_acceptor.spectrum(3:6); S_auto.spectrum(3:6)];
            
    %unmix =%%receive all non-fret concentrations
        t_concentrations=nonnegative_unmix(S,X,'verbose',true);
        

        
        
%% multiply all non-fret concentrations by correscponding full specrta and
%sum
    full_spectra=[S_donor.spectrum; S_acceptor.spectrum; S_auto.spectrum];
    t_reconstructed_bleach_spectrum=t_concentrations*full_spectra;
    
%% subtract the sum from the bleaching spectrum
    t_residuals_fret=S_bleach-t_reconstructed_bleach_spectrum;
    
%% divide this difference by acceptor and donor concentrations to obtain FRET estimate
    t_acc_multiply_donor=t_concentrations(:,1).*t_concentrations(:,2);
    t_fret_estimate=t_residuals_fret./repmat(t_acc_multiply_donor,1,size(t_residuals_fret,2));

%see where (in time) this difference is constant, and average those - this
%would be your FRET-/+ spectrum
%% SUPPORT
% [C,ia,ib] = intersect(A,B,'rows') 
% also returns index vectors ia and ib, 
% such that C = A(ia,:) and C = B(ib,:).



