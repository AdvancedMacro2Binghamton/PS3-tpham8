% PARAMETERS
beta = .9932; %discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)
y_s = [1, b]; % endowment in employment states
PI = [.97 .03; .5 .5]; % transition matrix


% ASSET VECTOR
a_lo = -2; %lower bound of grid points
a_hi = 5;%upper bound of grid points
num_a = 200;

a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR q
q_min = 0.98;
q_max = 1.1;


% ITERATE OVER ASSET PRICES
aggsav = 1 ;
while abs(aggsav) >= 0.01 ;
    q_guess = (q_min + q_max) / 2;
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, a', q_guess * a);
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons<0)=-Inf
    
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(2, num_a);
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol >.000001;
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        v=ret+beta*repmat(permute((PI*v_guess),[3 2 1]),[num_a 1 1]);
         % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vfn, pol_idx]=max(v,[],2);
        v_tol=max(max(abs(permute(vfn,[3 1 2])-v_guess)));
        v_guess=[vfn(:,:,1)';vfn(:,:,2)'];
    end;
    % KEEP DECSISION RULE
    pol_idx=permute(pol_idx,[3 1 2])
    pol_fn = a(pol_idx);
    
    
    Mu=ones(size(pol_fn))/numel(pol_fn);
    Mu_tol=1
    while Mu_tol>=1e-0.6;
        [emp_ind, a_ind, mass]=find(Mu);
        
        MuNew=zeros(size(Mu));
        
        for ii=1:length(emp_ind)
            apr_ind=pol_idx(emp_ind(ii),a_ind(ii));
            MuNew (:,apr_ind)=MuNew(:,apr_ind)+ (PI(emp_ind(ii),:)*mass(ii))';
        end
        Mu_tol=max(max(abs(MuNew-Mu)));
        Mu=MuNew;
    end
        aggsav=sum(sum(Mu.*pol_fn));
        if aggsav>0
            q_min=q_guess;
        else
            q_max=q_guess;
        end
end
            
figure(1);
subplot(2,1,1)
plot(a,vfn(:,:,1),a,vfn(:,:,2)),legend('Employed Vfn','Unemployed Vfn'),title('Value Function')
subplot(2,1,2)
plot(a,pol_fn(1,:),a,pol_fn(2,:)),legend ('Employed pol_fn','Unemployed pol_fn'),title('Policy Function')

%%%Lorenz Curve Using Gini Package in Matlab
population=[Mu(2,:);Mu(1,:)];

wealth=[bsxfun(@plus,a,y_s(2));bsxfun(@plus,a,y_s(1))];

wealth(wealth<0)=0;



earnings=[repmat(y_s(1),num_a,1);repmat(y_s(2),num_a,1)]



figure(2);
subplot(1,2,1);
wealth_gini=gini(population,wealth,true); title('Wealth Lorenz Curve')
subplot(1,2,2);
earnings_gini=gini(population,earnings,true); title('Earnings Lorenz Curve')





        
        
        
    