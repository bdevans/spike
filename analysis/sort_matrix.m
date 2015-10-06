function [Mat] = sort_matrix(Mat)

% Function to sort a matrix by columns
% IN: Mat - the matrix
% OUT: Mat - the sorted matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nrow=size(Mat,1);
% % [nrow,ncol]=size(Mat);
% colweights=(1:nrow).^2;
colweights=(1:size(Mat,1)).^2;
v_sum=colweights*Mat;
lMat=[v_sum;Mat];
% tmp=sortrows(lMat');
% tmp=tmp';
tmp=sortrows(lMat')';
Mat=tmp(2:end,:);

% v_ind=(1:ncol);
% V=[v_sum;v_ind];
% tmp=sortrows(V');
% V=tmp';
% 
% Mat(:,V(2,j))
% 
% for i=1:row
%     v_index=(1:row);
%     v_sum=zeros(1,row);
%     for j=1:col
% %         k=(i*row)+j;
% %         if Mat(k)~=-1
% %             v_sum(i)=v_sum(i)+(Mat(k)*j*j);
% %         end
%         v_sum(i)=v_sum(i)+(Mat(i,j)*j*j);
%     end
% end
% 
% for i=1:row
%     for j=i:col
%         if v_sum(i)>v_sum(j)
%             temp_sum=v_sum(i);
%             v_sum(i)=v_sum(j);
%             v_sum(j)=temp_sum;
% 
%             temp_index=v_index(i);
%             v_index(i)=v_index(j);
%             v_index(j)=temp_index;
%         end
%     end
% end
% 
% for j=1:col
%     for i=1:row
% %         v_sum(i)=Mat((v_index(i)*row)+j);
%         v_sum(i)=Mat(i,j);
%     end
%     for i=1:row
% %         Mat((i*row)+j)=v_sum(i);#
%         Mat(i,j)=v_sum(i);
%     end
% end