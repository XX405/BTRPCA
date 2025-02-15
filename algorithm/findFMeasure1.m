function [fmax, S, bestTH,pre, rec] = findFMeasure1(tenX2, S0)
%% S0为真实图片，tenX2为实际前景
    fmax = 0;
    Vth = linspace(min(abs(tenX2(:))), max(abs(tenX2(:))), 5000);
    for idx = 1:length(Vth)
        th = Vth(idx);
        Stmp = abs(tenX2) > th;
        [pre, rec] = countPR(S0, Stmp);
        ftmp = 2*rec*pre/(pre+rec);%%ftmp是计算f值
        if ~isnan(ftmp) && ftmp > fmax
            fmax = ftmp;
            S = Stmp;
            bestTH = th;
        end
    end
    th =bestTH ;
    Stmp = abs(tenX2) > th;
    [pre, rec] = countPR(S0, Stmp);
    
    if fmax == 0
        S = Stmp;
        bestTH = 0;
    end
end


% function         [pre, rec] = countPR(S0, Stmp)
% CS = (S0.*Stmp);
% pre = sum(CS(:))/sum(Stmp(:));
% rec = sum(CS(:))/sum(S0(:));
% end