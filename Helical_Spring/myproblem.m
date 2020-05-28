function [x fval history] = myproblem(x0)
    history = [];
       
%     options = optimset('Display','iter');
    options = optimset('OutputFcn', @myoutput,'Display','iter');    
    fun = @w_numerical_function;
    [x,fval,exitflag,output] = fminsearch(fun,x0,options);
 
    function stop = myoutput(x,optimvalues,state);
        stop = false;
        if isequal(state,'iter')
          history = [history; x];
        end
    end
    
%     function z = objfun(x)
%       z = exp(x(1))*(4*x(1)^2+2*x(2)^2+x(1)*x(2)+2*x(2));
%     end
end