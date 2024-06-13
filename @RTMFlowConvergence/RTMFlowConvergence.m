classdef RTMFlowConvergence < RTMFlow
    methods

        function [sample_obj, obj] = run(obj,sample_time)

          sample_obj = [];

          while ~obj.is_fully_saturated()


             %% Solve pressure & velocity Problem
             obj.pressure_class = obj.pressure_class.solve();
             obj.velocity_class = obj.velocity_class.compute_velocity(obj.pressure_class);
            
              
             %% Solve flow problem
             obj = obj.compute_flow_rates();
            
             %% Visualise
             obj.visualise_class.plot(obj);

             %% save simulation at first instance over time point
             if obj.pressure_class.time >= sample_time && isempty(sample_obj)
                 sample_obj = obj;
             end
            
             %% Increment to new time
             obj = obj.update_time_level();
         
             %% Update flow volumes and moving boundaries
             obj = obj.update_filling_percentage();
             obj = obj.update_computational_domain();

             
          end

          disp("end")

        end
    end
end