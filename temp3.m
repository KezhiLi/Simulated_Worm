len = zeros(1,size(all_skeletons,2));

for ii=1:size(all_skeletons,2)
    curr_ske = all_skeletons{ii};
  
    diff_curr_ske = curr_ske(2:end,:)-curr_ske(1:end-1,:);
    len(ii) = sum((diff_curr_ske(:,1).^2+diff_curr_ske(:,2).^2).^0.5);
end