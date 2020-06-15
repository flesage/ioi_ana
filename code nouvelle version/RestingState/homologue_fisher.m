function homologue_z = homologue_fisher(corr)
    fish = 0.5.*log((1+corr)./(1-corr));
    for i=1:12
        if (i<7)
            homologue(i,1) = fish(i,i+6);
        else
            homologue(i-6,2) = fish(i,i-6);
        end
    end
homologue_z = mean(homologue,2);
end

