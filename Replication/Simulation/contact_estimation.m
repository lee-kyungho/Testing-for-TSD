
function B = contact_estimation(r_N,D_collection,A,c_N)

% ---- D_collection --
% Dimension: n_grid * J * b
% ---- A -------------
% Dimension <= J - vecotor

N_J = 1:size(D_collection,2);
A_complement = setdiff(N_J,A); 

D_collection(:,A_complement) = 0;

% We pick index

Contact_A      = sum(abs(r_N.*D_collection(:,A)) <= c_N,2) == size(A,2);
Contact_A_comp = sum(r_N.*D_collection(:,A_complement) < -c_N,2) == size(A_complement,2);
Contact        = sum(Contact_A+Contact_A_comp,2) == 2;

% (r_N*D_collection(:,A_complement) <= c_N);

B = Contact;

end
