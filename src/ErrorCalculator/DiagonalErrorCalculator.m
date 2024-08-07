classdef DiagonalErrorCalculator < ErrorCalculator

    methods
        function obj = DiagonalErrorCalculator(data_curve, modeled_curve)
           obj = obj@ErrorCalculator(data_curve, modeled_curve);
        end

        function e = get_error(obj)
            obj.X = obj.X(1:end);
            obj.Y = obj.Y(1:end);
            Ydat = obj.Y(1:end-1);
            Xdat = obj.X(1:end-1);
            e = 0;
            for i = 1:length(Ydat)
                delta_x = abs(interp1(obj.Yhat, obj.Xhat, Ydat(i), 'linear', 'extrap') - Xdat(i))/max(obj.X);
                delta_y = abs(interp1(obj.X, obj.Y, interp1(obj.Yhat, obj.Xhat, Ydat(i), 'linear', 'extrap'), 'linear', 'extrap') - Ydat(i))/max(obj.Y);
                delta_o = delta_y*cos(atan(delta_y/delta_x));
                e = e + (delta_o^2);
            end
            e = sqrt(e)/length(obj.Y);
        end
    end
end