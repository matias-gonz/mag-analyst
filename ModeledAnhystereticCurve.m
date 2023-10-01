classdef ModeledAnhystereticCurve

    properties
        H
        M
        Mi
        dMdH
        dMidH
        HdMdH
        HdMidH
    end

    methods (Access = public)
        function obj = ModeledAnhystereticCurve(H, a, alpha, alphaMs, Ms)
            obj.H = H;
            [obj.M, obj.Mi] = obj.get_M(H, a, alphaMs, Ms);
            [obj.dMdH, obj.dMidH] = obj.get_dMdH(H, alpha, Ms, a, alphaMs);
            [obj.HdMdH, obj.HdMidH] = obj.get_HdMdH(H, dMidH);
        end
    end

    methods (Access = private)

        function [M, Mi] = get_M(H, a, alphaMs, Ms)
            M = zeros(1, length(H));
            Mi = zeros(length(a), length(H));
    
            for i = 1:length(H)
                for j = 1:length(a)
                    m = Utils().get_m(H(i), a(j), alphaMs(j));
                    Mi(j, i) = Ms(j)*m;
                    M(i) = M(i) + Mi(j, i);
                end 
            end
        end

        function [dMdH, dMidH] = get_dMdH(H, alpha, Ms, a, alphaMs)
            dMdH = zeros(1, length(H));
            dMidH = zeros(length(a), length(H));
        
            for i = 1:length(H)
                for j = 1:length(a)
                    m = Utils().get_m(H(i), a(j), alphaMs(j));
                    h = (H(i)+m*alphaMs(j))/a(j);
                    numerator = Ms(j)*Langevin(h,1)/a(j);
                    dMidH(j, i) = numerator/(1-alpha(j)*numerator);
                    dMdH(i) = dMdH(i) + dMidH(j, i);
                end  
            end
        end

        function [HdMdH, HdMidH] = get_HdMdH(H, dMidH)
            HdMdH = zeros(1, length(H));
            HdMidH = zeros(size(dMidH,1), length(H));
        
            for i = 1:length(H)
                for j = 1:size(dMidH,1)
                    HdMidH(j, i) = H(i)*dMidH(j, i);
                    HdMdH(i) = HdMdH(i) + HdMidH(j, i);
                end
            end
        end
    end
end