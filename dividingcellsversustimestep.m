
% Finding division and delamnation rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Just TACs
% resultTAC( isnan(resultTAC(:,1)), : ) = [];
% resultTAC( isinf(resultTAC(:,1)), : ) = [];
% cumaver = cumsum(resultTAC(:,1)) ./ (1:length(resultTAC(:,1)))';
% 
% plot(resultTAC(:,2),resultTAC(:,1))
% 
% hold on
% 
% plot(resultTAC(:,2),cumaver)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Just LESCs
% resultLESC( isnan(resultLESC(:,1)), : ) = [];
% resultLESC( isinf(resultLESC(:,1)), : ) = [];
% cumaver = cumsum(resultLESC(:,1)) ./ (1:length(resultLESC(:,1)))';
% 
% plot(resultLESC(:,2),resultLESC(:,1))
% 
% hold on
% 
% plot(resultLESC(:,2),cumaver)

% % 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Just Delamination
resultex( isnan(resultex(:,1)), : ) = [];
resultex( isinf(resultex(:,1)), : ) = [];
cumaver = cumsum(resultex(:,1)) ./ (1:length(resultex(:,1)))';

plot(resultex(:,2),resultex(:,1))

hold on

plot(resultex(:,2),cumaver)




