function [output,flag] = expm_1sing(A,dt,tol,flag)
% Function designed to compute A^(-1) (exp(Adt) -I) using Taylor series
% does this by computing the sum, until the next term below tolerance
% 'tol', then truncate the sum

% note flag = 1 means we had to use Taylor polynomial. flag = 0 means not
if isdiag(A)
    if flag == 0 %force direct method
        B = diag(A);
        nonz = find(B); %nonzero entries
        B_temp = diag(B(nonz)); %no 0's on diag, now invertible
        solved = B_temp\(expm(B_temp*dt) - eye(size(B_temp)));
        putin = ones(size(A,1),1); putin(nonz) = diag(solved);
        putin(B==0) = dt; %put back in the elements
        output = diag(putin);
    elseif flag == 1 %force Taylor poly method
        k = 1;
        sum = spzeros(size(A,1),1); %running sum
        curr = speye(size(A,1))*dt; %current term
        while norm(diag(curr)) > tol
            sum = sum + curr;
            k=k+1; curr = (A*dt/k)*curr; %update current term
        end
        output = sum;
    else %do preferred
        if max(max(A*dt)) > 13.5*pi %heuristic from observations
            disp('preferred is direct');
            B = diag(A);
            nonz = find(B); %nonzero entries
            B_temp = diag(B(nonz)); %no 0's on diag, now invertible
            solved = B_temp\(expm(B_temp*dt) - eye(size(B_temp)));
            putin = ones(size(A,1),1); putin(nonz) = diag(solved);
            putin(B==0) = dt; %put back in the elements
            output = diag(putin); flag = 0;
        else
            disp('preferred is Taylor');
            k = 1;
            sum = zeros(size(A,1),1); %running sum
            curr = speye(size(A,1))*dt; %current term
            while norm(diag(curr)) > tol
                sum = sum + curr;
                k=k+1; curr = (A*dt/k)*curr; %update current term
            end
            output = sum; flag = 1;
        end
    end
else %if not diagonal, there is no other method, just accept fate.
    k = 1;
    sum = zeros(size(A,1),1); %running sum
    curr = eye(size(A,1))*dt; %current term
    while norm(curr) > tol
        sum = sum + curr;
        k=k+1; curr = (A*dt/k)*curr; %update current term
    end
    output = sum; flag = 1;
end
end