function [Obj, B, Theta] = func_OPT_3diag(hRT_norm, hIT_norm, hRI_norm, NG)
% The received power is returned

NI = size(hIT_norm,1); % Number of RIS elements
if NG == 0 % Fully connected
    NG = NI;
end
G = NI/NG; % Number of groups


if NG == 1 % Single connected

    % Compute Theta
    theta = - angle(hRI_norm) - angle(hIT_norm.');
    Theta = diag(exp(1i * theta));

    % Compute B
    B = imag(2 * inv(Theta + eye(NI))) / 50;

else % Group or fully connected

    B = [];
    for g = 1:G
    
        % Truncated channels
        hRI_g = hRI_norm(NG*(g-1)+1:NG*g);
        hIT_g = hIT_norm(NG*(g-1)+1:NG*g);
        hRI_g_norm = hRI_g / norm(hRI_g);
        hIT_g_norm = hIT_g / norm(hIT_g);

        % Vectors alpha and beta
        a = 1i * 50 * (hRI_g_norm' * exp(1i * angle(hRT_norm)) + hIT_g_norm);
        b = hIT_g_norm - hRI_g_norm' * exp(1i * angle(hRT_norm));
        
        % Matrix A
        A = zeros(2*NG,3*NG);
        
        % First N columns
        A(1:NG,     1:NG)=diag(real(a));
        A(NG+1:2*NG,1:NG)=diag(imag(a));
        
        % Second N-1 columns
        A(1:NG,     NG+1:2*NG)=diag(real([a(2:NG);0])) + diag(real(a(1:NG-1)),-1);
        A(NG+1:2*NG,NG+1:2*NG)=diag(imag([a(2:NG);0])) + diag(imag(a(1:NG-1)),-1);
        
        A = A(:,1:2*NG-1);
        
        % Solve the linear system
        x = A \ [real(b);imag(b)]; % Same as pinv(A) * [real(b);imag(b)] and inv(A'*A)*A' * [real(b);imag(b)]
        
        % Matrix B
        D1 = x(1:NG);
        D2 = x(NG+1:end);
        B_tmp = diag(D1) + diag(D2,-1) + diag(D2,+1);
        
        B = blkdiag(B,B_tmp); % B is block diagonal

    end

    Theta = (eye(NI) + 1i*50*B) \ (eye(NI) - 1i*50*B); % Theta matrix
    
end

Obj = norm(hRT_norm + hRI_norm*Theta*hIT_norm) ^ 2; % Received Power

end