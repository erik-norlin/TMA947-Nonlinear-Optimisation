%------------------------------------------------------------
function fmincon1
  LB = [-inf; -inf];     % lower and upper bounds on the variables
  UB = [inf; inf];    % inf means infinity

  A = [-2,1]; b = [0];
  Aeq = []; beq = [ ]; % we don't have any linear inequality
                      % or equality constraints

  x0 = [0; 0];        % starting point

  % now, set some options of the solver
  % for more options, type 'help optimset' in Matlab
  options = optimset( 'LargeScale', 'off', ...
                      'Display', 'iter', ...
	              'GradObj', 'on', ...
                      'GradCon', 'on', ...
                      'DerivativeCheck', 'on' );

  
  % run the solver
  [x,z,ef,output,lambda] = fmincon( @testf, x0, ...
                                    A, b, ...
                                    Aeq, beq, ...
                                    LB, UB, ...
                                    @testg, options );
  
  % if success, print out the solution and multipliers
  if( ef > 0 )
    disp( sprintf( '============================================================' ) );
    disp( sprintf( 'Optimal solution: ') ); disp( x' );
    disp( sprintf( 'Optimal value:    %g', z ) );
    disp( 'Lagrange multipliers: ' );
    disp( sprintf( 'Upper bounds    :' ) ); disp( lambda.lower' );
    disp( sprintf( 'Lower bounds    :' ) ); disp( lambda.upper' );
    disp( sprintf( 'Lin. ineq.      :' ) ); disp( lambda.ineqlin' );
    disp( sprintf( 'Lin. eq.        :' ) ); disp( lambda.eqlin' );
    disp( sprintf( 'Nonlin. ineq.   :' ) ); disp( lambda.ineqnonlin' );
    disp( sprintf( 'Nonlin. eq.     :' ) ); disp( lambda.eqnonlin' );
  end;

  return;

%------------------------------------------------------------
function [f, df] = testf(x)
  % Calculate the objective function value f and
  % its gradient df

  f = -x(1) + 2*x(1)^2 - 2*x(2) + x(2)^2 - x(1)*x(2);
  df = [-1+4*x(1)-x(2); -2+2*x(2)-x(1)];
  return;

%------------------------------------------------------------
function [C, Ceq, dC, dCeq] = testg(x)
  % Calculate the values of the nonlinear inequality
  % constraints C(x) <= 0, equality constraints Ceq(x) = 0
  % and the corresponding gradients.

   C    = [x(1)^2-x(2)]; 
   Ceq  = []; % no equality constraints
   dC   = [[2*x(1); -1]];
   dCeq = [];

   return;
