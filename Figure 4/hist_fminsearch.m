function [x, fval, history, obj_val_his] = hist_fminsearch(my_wanted_params, LB_single, objfunc_single, max_fcn, max_iter)
% Save parameter values and cost for each iteration in fminsearchbnd
history     = [];
obj_val_his = [];
options = optimset('OutputFcn', @myoutput, 'MaxFunEvals', max_fcn, 'MaxIter', max_iter);

[x, fval] = fminsearchbnd(objfunc_single, my_wanted_params,LB_single,[],options);

    function stop = myoutput(x,optimvalues,state)
        stop = false;
        if isequal(state,'iter')
            history = [history; x];
            obj_val_his = [obj_val_his; optimvalues.fval];
            itr_disp = ['Iteration: ' num2str(optimvalues.iteration)];
            disp(itr_disp)
        end
    end

end