classdef performancemlapp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        LoadCovariateDataButton      matlab.ui.control.Button
        NumberofKnotsEditField       matlab.ui.control.NumericEditField
        NumberofKnotsEditFieldLabel  matlab.ui.control.Label
        SaveResultsButton            matlab.ui.control.Button
        PlotResultsButton            matlab.ui.control.Button
        ThinningEditField            matlab.ui.control.NumericEditField
        ThinningEditFieldLabel       matlab.ui.control.Label
        NumberofIterationsEditField  matlab.ui.control.NumericEditField
        NumberofIterationsEditFieldLabel  matlab.ui.control.Label
        BurnInEditField              matlab.ui.control.NumericEditField
        BurnInEditFieldLabel         matlab.ui.control.Label
        RunCodeButton                matlab.ui.control.Button
        LoadPerformanceDataButton    matlab.ui.control.Button
    end

    
    properties (Access = private)
        Property % Description
        numberofits = 100;
        numberofknots = 100;
        burnin = 100;
        every = 1;
        fig1 = uifigure('Name','Dialog Box','Position',[100 100 640 175]);
        PathName = [];
        datamat_pred = [];
        ind1_pred = [];
        data_athlete = [];
        data_performance = [];
        prob_inclusion = [];
        alpha = [];
        gamma = [];
        beta = [];
        beta0 = [];
        sigmasq = [];
        output = [];
        sigmasqalpha = [];
        sigmasqbeta = [];
        knots = [];
        age_knots = [];
        age_xstar = [];
        IDstar = [];
        fig2 = uifigure('Name','Plot Box','Position',[100 100 800 800]);
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: LoadPerformanceDataButton
        function LoadPerformanceDataButtonPushed(app, event)
            [filename, pathname] = uigetfile({'*.csv'},'FileSelector');
            if ( ~isequal(filename, 0) )
                app.PathName = pathname;
                name1 = [pathname filename];
                app.data_performance = readtable(name1);
                
                uialert(app.fig1,'Performance Data Successfully Loaded','Load Data','Icon','success')
            end
        end

        % Value changed function: BurnInEditField
        function BurnInEditFieldValueChanged(app, event)
            value = app.BurnInEditField.Value;
            app.burnin = value;
        end

        % Value changed function: NumberofIterationsEditField
        function NumberofIterationsEditFieldValueChanged(app, event)
            value = app.NumberofIterationsEditField.Value;
            app.numberofits = value;
        end

        % Value changed function: ThinningEditField
        function ThinningEditFieldValueChanged(app, event)
            value = app.ThinningEditField.Value;
            app.every = value;
        end

        % Button pushed function: RunCodeButton
        function RunCodeButtonPushed(app, event)

            uialert(app.fig1,'Run Started','Run','Icon','success')
            
            n = size(app.data_performance, 1);
            
            data_fixed1 = app.data_athlete;
            data_local = zeros(n, 2);
            data_local(:, 1) = app.data_performance.Age;
            data_local(:, 2) = app.data_performance.Performance;
           
            
            datenumhold = data_local(:, 1);
            datenum1 = data_local(:, 1);
            datenum1 = (datenum1 - min(datenum1)) / (max(datenum1) - min(datenum1));
            
            app.knots = linspace(0, 1, app.numberofknots + 1);
            app.knots = app.knots(1:(end-1));
            app.age_knots = min(datenumhold) + (max(datenumhold) - min(datenumhold)) * app.knots;
            
            app.IDstar = unique(app.data_performance.ID);
            
            xstar = linspace(0, 1, 101);
            app.age_xstar = min(datenumhold) + (max(datenumhold) - min(datenumhold)) * xstar;
            data_all = cell(length(app.IDstar), 1);
            target_all = cell(length(app.IDstar), 1);
            data_fixed = cell(length(app.IDstar), 1);
            
            [app.datamat_pred, app.ind1_pred] = create_datamat(xstar, app.knots, 1);
                      
            for i = 1:length(app.IDstar)
                
                ID = app.IDstar(i);
                
                z = app.data_performance.ID==ID;
                
%                 sub2 = (app.age_xstar > min(data_local(z==1, 1))) .* (app.age_xstar < max(data_local(z==1, 1)));
%                 x1 = find(sub2==1, 1, 'first');
%                 sub2(x1 - 1) = 1;
%                 x1 = find(sub2==1, 1, 'last');
%                 sub2(x1 + 1) = 1;
                
                data = datenum1(z==1);
                %y = athleteSubset.MovingAverage(z==1);
                target = app.data_performance.Performance(z==1);
                
                [datamat, ind1] = create_datamat(data, app.knots, 1);
                
                data_all{i, 1} = datamat;
                target_all{i, 1} = target;
                data_fixed{i, 1} = data_fixed1(z, :);
            end
            
            [app.prob_inclusion, app.alpha, app.beta, app.gamma, app.beta0, app.sigmasq, app.output] = ASI_sampler_special(data_fixed, data_all, target_all, 100, 2, 2, 0, 0, [1 (app.numberofknots - 3) / 3], 0.234, 1, app.burnin, app.numberofits, app.every, 25, 1, 1, 1, app.fig1);
            
            
            uialert(app.fig1,'Run Completed','Run','Icon','success')

        end

        % Button pushed function: LoadCovariateDataButton
        function LoadCovariateDataButtonPushed(app, event)
            [filename, pathname] = uigetfile({'*.csv'},'FileSelector');
            if ( ~isequal(filename, 0) )
                app.PathName = pathname;
                name1 = [pathname filename];
                app.data_athlete = readmatrix(name1);
                
 
                
                uialert(app.fig1,'Athlete Data Successfully Loaded','Load Data','Icon','success')
            end
        end

        % Button pushed function: SaveResultsButton
        function SaveResultsButtonPushed(app, event)
            
          prob_inclusion = app.prob_inclusion;
          alpha = app.alpha;
          b = app.beta;
          included = app.gamma;
          gamma = app.beta0;
          sigmasqalpha = app.output.sigmasqalpha;
          sigmasqbeta = app.output.sigmasqbeta;
          knots = app.knots;
          sigmasq = app.sigmasq;
          age_knots = app.age_knots;
          age_xstar = app.age_xstar;
          datamat_pred = app.datamat_pred;
          ind1_pred = app.ind1_pred;
          
          [file,path] = uiputfile('*.mat','Workspace File');
          if ( ~isequal(file, 0) )
              name1 = [path file];
              
              save(name1,'prob_inclusion', 'included', 'alpha', 'b', 'gamma', 'sigmasqalpha', 'sigmasqbeta', 'knots', 'age_knots', 'age_xstar', 'datamat_pred', 'ind1_pred', 'sigmasq');
              
              if ( ishandle(app.fig1) == 0 )
                  app.fig1 = uifigure('Name','Dialog Box','Position',[100 100 640 175]);
              end
              
              uialert(app.fig1,'Results Successfully Saved','Save Results','Icon','success')
          end
        end

        % Button pushed function: PlotResultsButton
        function PlotResultsButtonPushed(app, event)

            if ( ishandle(app.fig2) == 0 )
                app.fig2 = uifigure('Name','Plot Box','Position',[100 100 800 800]);
            end
            ax_graph = uiaxes('Parent', app.fig2, 'Position', [100 100 600 600]);

            if ( isnumeric(app.IDstar))
                local_IDstar = cell(length(app.IDstar), 1);
                for i = 1:length(local_IDstar)
                    local_IDstar{i, 1} = num2str(app.IDstar(i));
                end
            else
                local_IDstar = app.IDstar; 
            end
            
            uilabel(app.fig2,'Text','Athlete ID', 'Position', [530 730 100 22])           
            dd_ID = uidropdown(app.fig2, 'Position',[600 730 100 22], 'Items', local_IDstar, 'ValueChangedFcn', @(dd_ID, event) update_graph(dd_ID, ax_graph));
            
            
            function update_graph(dd_orders, ax_graph)
                
                if ( ishandle(app.fig2) == 0 )
                    app.fig2 = uifigure('Name','Plot Box','Position',[100 100 800 800]);
                end
                ax_graph = uiaxes('Parent', app.fig2, 'Position', [100 100 600 600]);

                
                pos = find(strcmp(local_IDstar, dd_ID.Value));
                local_alpha = app.alpha(pos, :);
                local_beta = zeros(size(app.datamat_pred, 2), size(app.beta, 2));
                for i1 = 1:size(app.beta, 2)
                    local_beta(app.gamma{pos, i1}, i1) = app.beta{pos, i1};
                end
                predmean = ones(size(app.datamat_pred, 1), 1) * local_alpha + app.datamat_pred * local_beta;
                
                x = quantile(predmean, [0.025 0.975], 2);
                
                z = app.data_performance.ID == app.IDstar(pos);
                age_min = min(app.data_performance.Age(z == 1));
                age_max = max(app.data_performance.Age(z == 1));
                
                mask = (app.age_xstar > age_min) .* (app.age_xstar < age_max);
                
                plot(ax_graph, app.age_xstar(mask == 1), median(predmean(mask == 1, :), 2), 'b-')
                hold(ax_graph, 'on')
                plot(ax_graph, app.age_xstar(mask == 1), x(mask == 1, 1), 'b--')
                plot(ax_graph, app.age_xstar(mask == 1), x(mask == 1, 2), 'b--')
                hold(ax_graph, 'off')
                
                xlabel(ax_graph, 'Age')
                ylabel(ax_graph, 'Mean Standardized Performance')
                
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create LoadPerformanceDataButton
            app.LoadPerformanceDataButton = uibutton(app.UIFigure, 'push');
            app.LoadPerformanceDataButton.ButtonPushedFcn = createCallbackFcn(app, @LoadPerformanceDataButtonPushed, true);
            app.LoadPerformanceDataButton.Position = [95 387 143 22];
            app.LoadPerformanceDataButton.Text = 'Load Performance Data';

            % Create RunCodeButton
            app.RunCodeButton = uibutton(app.UIFigure, 'push');
            app.RunCodeButton.ButtonPushedFcn = createCallbackFcn(app, @RunCodeButtonPushed, true);
            app.RunCodeButton.Position = [47 181 100 22];
            app.RunCodeButton.Text = 'Run Code';

            % Create BurnInEditFieldLabel
            app.BurnInEditFieldLabel = uilabel(app.UIFigure);
            app.BurnInEditFieldLabel.HorizontalAlignment = 'right';
            app.BurnInEditFieldLabel.Position = [412 320 46 22];
            app.BurnInEditFieldLabel.Text = 'Burn-In';

            % Create BurnInEditField
            app.BurnInEditField = uieditfield(app.UIFigure, 'numeric');
            app.BurnInEditField.ValueChangedFcn = createCallbackFcn(app, @BurnInEditFieldValueChanged, true);
            app.BurnInEditField.Position = [478 320 72 22];
            app.BurnInEditField.Value = 100;

            % Create NumberofIterationsEditFieldLabel
            app.NumberofIterationsEditFieldLabel = uilabel(app.UIFigure);
            app.NumberofIterationsEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofIterationsEditFieldLabel.Position = [342 284 116 22];
            app.NumberofIterationsEditFieldLabel.Text = 'Number of Iterations';

            % Create NumberofIterationsEditField
            app.NumberofIterationsEditField = uieditfield(app.UIFigure, 'numeric');
            app.NumberofIterationsEditField.ValueChangedFcn = createCallbackFcn(app, @NumberofIterationsEditFieldValueChanged, true);
            app.NumberofIterationsEditField.Position = [478 284 72 22];
            app.NumberofIterationsEditField.Value = 100;

            % Create ThinningEditFieldLabel
            app.ThinningEditFieldLabel = uilabel(app.UIFigure);
            app.ThinningEditFieldLabel.HorizontalAlignment = 'right';
            app.ThinningEditFieldLabel.Position = [384 243 74 22];
            app.ThinningEditFieldLabel.Text = 'Thinning';

            % Create ThinningEditField
            app.ThinningEditField = uieditfield(app.UIFigure, 'numeric');
            app.ThinningEditField.ValueChangedFcn = createCallbackFcn(app, @ThinningEditFieldValueChanged, true);
            app.ThinningEditField.Position = [478 243 72 22];
            app.ThinningEditField.Value = 1;

            % Create PlotResultsButton
            app.PlotResultsButton = uibutton(app.UIFigure, 'push');
            app.PlotResultsButton.ButtonPushedFcn = createCallbackFcn(app, @PlotResultsButtonPushed, true);
            app.PlotResultsButton.Position = [271 181 100 22];
            app.PlotResultsButton.Text = 'Plot Results';

            % Create SaveResultsButton
            app.SaveResultsButton = uibutton(app.UIFigure, 'push');
            app.SaveResultsButton.ButtonPushedFcn = createCallbackFcn(app, @SaveResultsButtonPushed, true);
            app.SaveResultsButton.Position = [478 181 100 22];
            app.SaveResultsButton.Text = 'Save Results';

            % Create NumberofKnotsEditFieldLabel
            app.NumberofKnotsEditFieldLabel = uilabel(app.UIFigure);
            app.NumberofKnotsEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofKnotsEditFieldLabel.Position = [17 320 97 22];
            app.NumberofKnotsEditFieldLabel.Text = 'Number of Knots';

            % Create NumberofKnotsEditField
            app.NumberofKnotsEditField = uieditfield(app.UIFigure, 'numeric');
            app.NumberofKnotsEditField.Position = [138 320 56 22];
            app.NumberofKnotsEditField.Value = 100;

            % Create LoadCovariateDataButton
            app.LoadCovariateDataButton = uibutton(app.UIFigure, 'push');
            app.LoadCovariateDataButton.ButtonPushedFcn = createCallbackFcn(app, @LoadCovariateDataButtonPushed, true);
            app.LoadCovariateDataButton.Position = [405 387 126 22];
            app.LoadCovariateDataButton.Text = 'Load Covariate Data';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = performance

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end