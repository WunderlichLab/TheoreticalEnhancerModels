
function PatternPlotter(figTracker, yPointsOrig, yPointsOrigSpecial ...
                     , myColor, width, fixedString, fixedCount, enhancerBool)
         figure(figTracker)
         if enhancerBool
            yPointsOrig = reshape(yPointsOrig,[1,size(yPointsOrig,3)]);
         end
         if enhancerBool
            yPointsOrigSpecial = reshape(yPointsOrigSpecial,[1,size(yPointsOrigSpecial,3)]);
         end

         % finding non-zero indices which correspond
         % to the total number of binding sites for T1 or T2.
         % these are naturally the x-axis 
         xPoints = find(yPointsOrig)';
         % && ~isempty(yPointsOrig)
         if ~isempty(xPoints) && length(xPoints) > 1 
           % setting the y values to be the values at these points
           yPoints = yPointsOrig(xPoints)';
           % we do this after since at the entry i of the matrix
           % actually corresponds to i - 1 binding sites  
           if ~enhancerBool
            xPoints = xPoints - 1;                
           end
           plot(xPoints,yPoints,'LineWidth',width,'Color',myColor)
           hold on
           title(strcat(fixedString,' num: ',num2str(fixedCount)))

end