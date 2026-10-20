function [stored_metamodels,single_errors,time_data] = adaptive_sampling_process(metamodel_ini,M,adaptive_method,A_sampling,number_of_adaptive_iterations, opti_strategy, iteration_number)
% Details: Process of adding new samples until threshold reached
%
% inputs:
% metamodel_ini - initial metamodel
% M - respone function
% adaptive_method - adaptive sampling technique
% A_sampling - parameter space
% number_of_adaptive_iterations - maximum number of iterations
%
% outputs:
% stored_metamodels - stored metamodels
% single_errors - errors of each adaptive step


global plotFlag
global isFEAFlag
if isempty(isFEAFlag)
    isFEAFlag = false;
end
iter = 1;
stored_metamodels{1} = metamodel_ini;
metamodel = metamodel_ini;


lb = A_sampling(1,:);
ub = A_sampling(2,:);

lb_unity = zeros(size(lb));
ub_unity = ones(size(ub));
A_sampling_unity =[lb_unity;ub_unity];


Y = metamodel_ini.Y;
X = metamodel_ini.X;

    n_of_Variables = size(metamodel_ini.X,2);
    no_test_points = 5000 * n_of_Variables;
    test_points_unscaled = lhs_scaled(no_test_points,A_sampling(1,:),A_sampling(2,:));
    if n_of_Variables == 1
    test_points_unscaled = sort(test_points_unscaled);
    end
    test_points = scale_vector_to_unity(lb, ub, test_points_unscaled);
if plotFlag == 1


    for i=1:no_test_points
        test_points_response(i) = M(test_points_unscaled(i,:));
    end
    
    minTestPoints = min(test_points_response);
    maxTestPoints = max(test_points_response);
    
    maxTestPointsLimit = maxTestPoints + abs(maxTestPoints)/10;
    minTestPointsLimit = minTestPoints - abs(minTestPoints)/10;
    minTestPointsLimitUpper = minTestPoints - abs(minTestPoints)/20;
    
    
    
    
    for i=1:no_test_points
        [metamodel_response(i),~] = metamodel_ini.predict(test_points(i,:));
    end
    
    if n_of_Variables == 1
        figure(1)
        plot(test_points, test_points_response,'--','LineWidth',2.0); hold on;
        plot(test_points,metamodel_response ,'LineWidth',2.0); hold on;
        scatter(X, minTestPointsLimitUpper*ones(numel(X),1),30,'k', 'filled'); hold on;
        hold off;
        %               scatter(x_new, -10,60,'k', 'filled'); hold on;
        
        xlabel('$x$', 'Interpreter', 'Latex','FontSize',18);
        ylabel('$y$', 'Interpreter', 'Latex','FontSize',18);
        ST1 = ['$M$'];
        ST2 = ['$\hat{M}$'];
        ldg = legend({ ST1, ST2},'Interpreter', 'Latex','FontSize',12,'Location','north');
        ylim([minTestPointsLimit, maxTestPointsLimit])
        set(gca,'FontSize',18)
        drawnow()
    elseif n_of_Variables == 2
        xd = (max(test_points(:,1))-min(test_points(:,1)))/200;
        yd = (max(test_points(:,2))-min(test_points(:,2)))/200;
        [xq,yq] = meshgrid(min(test_points(:,1)):xd:max(test_points(:,1)), min(test_points(:,2)):yd:max(test_points(:,2)));
        
        
        h = figure(1) ;
        set(gcf, 'Renderer', 'painters', 'Position',[0 200 700 500]);
        subplot(2,1,1);
        [xq,yq] = meshgrid(min(test_points(:,1)):xd:max(test_points(:,1)), min(test_points(:,2)):yd:max(test_points(:,2)));
        vq = griddata(test_points(:,1),test_points(:,2),metamodel_response,xq,yq);

        s = surf(xq,yq,vq); hold on;
        hold off;
        xlabel('$x$', 'Interpreter', 'Latex','FontSize',18);
        ylabel('$y$', 'Interpreter', 'Latex','FontSize',18);
        zlabel('$\hat{M}$', 'Interpreter', 'Latex','FontSize',18);
        s.EdgeColor = 'none';
        
        ylim([0, 1])
        zlim([minTestPointsLimit, maxTestPointsLimit])
        xlim([0,1])
        set(gca,'FontSize',18)
        
        error =   abs(metamodel_response-    test_points_response)/(abs(max(test_points_response)-min(test_points_response)));
        min_error = 0;
        max_error = max(error);

        hold off;
        
        subplot(2,1,2);
        
        vq = griddata(test_points(:,1),test_points(:,2),error,xq,yq);
        s = surf(xq,yq,vq); hold on;
        scatter3(X(:,1),X(:,2), 1000*ones(size(X(:,1))),80,'k', 'filled'); hold on;
        
        xlabel('$x$', 'Interpreter', 'Latex','FontSize',18);
        ylabel('$y$', 'Interpreter', 'Latex','FontSize',18);
        %     zlabel('$M_{Bu,2d}$', 'Interpreter', 'Latex','FontSize',18);
        s.EdgeColor = 'none';
        s.FaceAlpha = 0.9;
        ylim([0, 1])
        xlim([0,1])
        colormap(jet)
        colorbar
        caxis([min_error ,max_error ])
        %       h( 'Position', [10 10 900 600])
        %       h.Position = [10 10 700 600];
        set(gcf, 'Renderer', 'painters', 'Position',[10 10 700 600]);
        view(2)
        set(gca,'FontSize',18)
        
        hold off;
        drawnow()
    end
    
end

data_errors = Error_Saver();


if strcmp(adaptive_method,'AME')
    global gamma_index
    global AMEpattern
    AMEpattern = [0,0.5,1.0,100];
    gamma_index = 1;
elseif strcmp(adaptive_method,'MEPE')
    global MEPE_q
    MEPE_q = 1;
elseif strcmp(adaptive_method,'MTE')
    global MEPE_q
    MEPE_q = 1;
end

time_per_iter = zeros(number_of_adaptive_iterations,1);      % 每轮总耗时
time_get_sample = zeros(number_of_adaptive_iterations,1);    % 获取新采样点耗时
sample_number = zeros(number_of_adaptive_iterations,1);      % 当前样本数量


while (iter <= number_of_adaptive_iterations)

    t_iter = tic;   % <-- 开始计时
    clear y_new x_new

        % --------- (A) LHS_rebuild：每轮重建训练集 ----------
    if strcmp(adaptive_method,'LHS_rebuild')
        m0 = size(metamodel_ini.X,1);
        m_train = m0 + iter;

        % 1) 重建训练集：新的 LHS（物理空间）
        t_sample = tic;
        
        train_points_unscaled = lhs_scaled(m_train, A_sampling(1,:), A_sampling(2,:));
        
        time_get_sample(iter) = toc(t_sample);

        % 2) 缩放到 unity 给 OK_model 用
        X = scale_vector_to_unity(A_sampling(1,:), A_sampling(2,:), train_points_unscaled);

        % 3) 真实响应
        Y = zeros(m_train,1);
        for ss=1:m_train
            Y(ss,1) = M(train_points_unscaled(ss,:));
        end

        % 4) 训练新代理
        metamodel = OK_model(metamodel.auto_correlation_function, X, Y, metamodel.theta_opti_technique);

    else
        t_sample = tic;   % 开始记录“获取新采样点”的时间
        
        % 原逻辑：每次新增一个点并重训
        x_new = metamodel.adaptive_sampling(adaptive_method,A_sampling_unity, opti_strategy, iter, number_of_adaptive_iterations);

        time_get_sample(iter) = toc(t_sample);   % 结束记录
        
        x_new_scaled = scale_vector_from_unity(lb, ub, x_new);

        for ss=1:size(x_new_scaled,1)
            y_new(ss,1) = M(x_new_scaled(ss,:));
        end

        Y = [Y; y_new];
        X = [X; x_new];

        metamodel = OK_model(metamodel.auto_correlation_function, X, Y, metamodel.theta_opti_technique);
    end
    
    
    %iter = iter+1;
    
    %% Errors
    if ~isFEAFlag
        for i=1:no_test_points
            test_points_response(i) = M(test_points_unscaled(i,:));
        end

        for i=1:no_test_points
            [metamodel_response(i),~] = metamodel.predict(test_points(i,:));
        end
    end
    
    
    if plotFlag == 1
        if n_of_Variables == 1
            figure(1)
            plot(test_points, test_points_response,'--','LineWidth',2.0); hold on;
            plot(test_points,metamodel_response ,'LineWidth',2.0); hold on;
            scatter(X, minTestPointsLimitUpper*ones(numel(X),1),30,'k', 'filled'); hold on;
            scatter(x_new, minTestPointsLimitUpper,60,'r', 'filled'); hold off;
 xlabel('$x$', 'Interpreter', 'Latex','FontSize',18);
        ylabel('$y$', 'Interpreter', 'Latex','FontSize',18);
        ST1 = ['$M$'];
        ST2 = ['$\hat{M}$'];
        ldg = legend({ ST1, ST2},'Interpreter', 'Latex','FontSize',12,'Location','north');
        ylim([minTestPointsLimit, maxTestPointsLimit])
        set(gca,'FontSize',18)
            drawnow()
        elseif n_of_Variables == 2
            vq1 = griddata(test_points(:,1),test_points(:,2),metamodel_response',xq,yq);
            h = figure(1) ;
            
            subplot(2,1,1);
            % cla(h1);
            
            s = surf(xq,yq,vq1); hold on;
            s.EdgeColor = 'none';
            xlabel('$x$', 'Interpreter', 'Latex','FontSize',18);
        ylabel('$y$', 'Interpreter', 'Latex','FontSize',18);
        zlabel('$\hat{M}$', 'Interpreter', 'Latex','FontSize',18);
        ylim([0, 1])
        zlim([minTestPointsLimit, maxTestPointsLimit])
        xlim([0,1])
        set(gca,'FontSize',18)
            hold off;
            
            
            error =   abs(metamodel_response-    test_points_response)/(abs(max(test_points_response)-min(test_points_response))) ;
            
            vq2 = griddata(test_points(:,1),test_points(:,2),error',xq,yq);
            subplot(2,1,2);
            
            
            s = surf(xq,yq,vq2); hold on;
            scatter3(X(:,1),X(:,2), 1000*ones(size(X(:,1))),80,'k', 'filled'); hold on;
            scatter3(x_new(:,1),x_new(:,2), 1200,150,'r', 'filled'); hold on;
  xlabel('$x$', 'Interpreter', 'Latex','FontSize',18);
        ylabel('$y$', 'Interpreter', 'Latex','FontSize',18);
        %     zlabel('$M_{Bu,2d}$', 'Interpreter', 'Latex','FontSize',18);
        s.EdgeColor = 'none';
        s.FaceAlpha = 0.9;
        ylim([0, 1])
        xlim([0,1])
        %colormap(jet)
        colorbar
        caxis([min_error ,max_error ])
        %       h( 'Position', [10 10 900 600])
        %       h.Position = [10 10 700 600];
        set(gcf, 'Renderer', 'painters', 'Position',[10 10 700 600]);
        view(2)
        set(gca,'FontSize',18)
            hold off;
            
            drawnow()
        end
    end
    
    
    if ~isFEAFlag
        data_errors=data_errors.update(metamodel.m,test_points_response , metamodel_response );

        ST = [adaptive_method,'   Iteration_no: ', num2str(iteration_number), '   m = ', num2str(size(X,1)),'   NMAE:  ', num2str(data_errors.error_data.MAE(end))];
        disp(ST);
        ST =[adaptive_method,'   Iteration_no: ', num2str(iteration_number), '   m = ', num2str(size(X,1)),'   NMaxAE:  ', num2str(data_errors.error_data.RMAE(end))];
        disp(ST);
        ST =[adaptive_method,'   Iteration_no: ', num2str(iteration_number), '   m = ', num2str(size(X,1)),'   NRMSE:  ', num2str(data_errors.error_data.RMSE(end))];
        disp(ST);
        ST =[adaptive_method,'   Iteration_no: ', num2str(iteration_number), '   m = ', num2str(size(X,1)),'   R^2:  ', num2str(data_errors.error_data.R_sq(end))];
        disp(ST);
    else
        ST = [adaptive_method,'   Iteration_no: ', num2str(iteration_number), '   m = ', num2str(size(X,1)),'   FE validation skipped'];
        disp(ST);
    end
    fprintf('\n');
    
    stored_metamodels{iter} = metamodel;
    single_errors = {adaptive_method,data_errors};

    time_per_iter(iter) = toc(t_iter);   % 本轮总耗时
    sample_number(iter) = size(X,1);     % 当前样本数量

    fprintf('[%s] Repeat %d, Adaptive iter %d: sample time = %.4f s, total time = %.4f s\n', ...
        adaptive_method, iteration_number, iter, ...
        time_get_sample(iter), time_per_iter(iter));

    iter = iter+1;
    
end

    adaptive_iter = (1:number_of_adaptive_iterations)';
    
    time_data = table( ...
        adaptive_iter, ...
        sample_number, ...
        time_get_sample, ...
        time_per_iter, ...
        'VariableNames', { ...
            'AdaptiveIteration', ...
            'SampleNumber', ...
            'TimeGetSample_s', ...
            'TotalIterationTime_s'} ...
    );

end
