function [output,flag] = expm_2sing(A,dt,tol,flag)
% Function designed to compute A^(-2) (exp(Adt) -I - Adt) using Taylor series
% does this by computing the sum, until the next term below tolerance
% 'tol', then truncate the sum

% flag indicates which methods was tried. flag = 0 means direct, flag=1
% means Taylor poly, flag = 2 means try preferred.
if isdiag(A)
    if flag == 0 %force direct method
        B = diag(A);
        nonz = find(B); %nonzero entries
        B_temp = diag(B(nonz)); %no 0's on diag, now invertible
        solved = (B_temp^2)\((expm(B_temp*dt) - eye(size(B_temp)) - B_temp*dt)/dt);
        putin = ones(size(A,1),1); putin(nonz) = diag(solved);
        putin(B==0) = dt/2; %put back in the elements
        output = diag(putin);
    elseif flag == 1 %force Taylor poly method
        k = 2;
        sum = zeros(size(A,1),1); %running sum
        curr = eye(size(A))*dt / k; %current term
        while norm(diag(curr)) > tol
            sum = sum + curr;
            k=k+1; curr = (A*dt/k)*curr; %update current term
        end
        output = sum;   
    else %do preferred
        if max(max(A*dt)) > 13.5*pi %heuristic from observations
            B = diag(A);
            nonz = find(B); %nonzero entries
            B_temp = diag(B(nonz)); %no 0's on diag, now invertible
            solved = (B_temp^2)\((expm(B_temp*dt) - eye(size(B_temp)) - B_temp*dt)/dt);
            putin = ones(size(A,1),1); putin(nonz) = diag(solved);
            putin(B==0) = dt/2; %put back in the elements
            output = diag(putin); flag = 0; 
        else
            k = 2;
            sum = zeros(size(A,1),1); %running sum
            curr = eye(size(A))*dt / k; %current term
            while norm(diag(curr)) > tol
                sum = sum + curr;
                k=k+1; curr = (A*dt/k)*curr; %update current term
            end
            output = sum; flag = 1;
        end
    end
else %if not diagonal, there is no other method, just accept fate.
    disp(norm(A*dt));
    k = 2;
    sum = zeros(size(A,1),1); %running sum
    curr = eye(size(A))*dt / k; %current term
    while norm(curr) > tol
        sum = sum + curr;
        k=k+1; curr = (A*dt/k)*curr; %update current term
    end
    output = sum; flag = 1;
end
end
