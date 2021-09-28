function create_pca_analysis_figures(DFF)
%CREATE_PCA_ANALYSIS_FIGURES Creates a couple of figures containing
%analysis via PCA of the DFF of components

rescaled = mat2gray(DFF,[0 .25]);
[coeff,~,~,~,explained]=pca(rescaled');

figure();
hold on;

pcs = coeff(:,1:3);

[~,max_pc] = max(abs(pcs),[],2);
idxs = [];
num_per_group = [];
for pc_idx = 1:3
    pos = find(max_pc==pc_idx & pcs(:, pc_idx)>=0);
    neg = find(max_pc==pc_idx & pcs(:, pc_idx)<0);
    [~,reordered_pos]=sortrows(pcs(pos,:),pc_idx,'descend');
    [~,reordered_neg]=sortrows(pcs(neg,:),pc_idx,'ascend');
    idxs=[idxs; pos(reordered_pos); neg(reordered_neg)];
    num_per_group=[num_per_group; length(pos); length(neg)];
end
line_pos_y = cumsum(num_per_group);
line_pos_y = line_pos_y(1:end-1);
idxs = flip(idxs);
line_pos_y=(length(idxs)-line_pos_y)+0.5;

%% subplot 1
ax1 = subplot(2,2,1);
imshow(DFF(flip(idxs),:));
set(ax1,'CLim',[0, 0.25]);
colormap('jet');
axis normal;
set(ax1,'XTick',[]);

%% subplot 2
ax2 = subplot(2,2,2);
colors = [162 40 30;
    220 128 129;
    60 72 147;
    190 186 218;
    61 143 66;
    135 181 67]/255.0;
barWidth = 0.75;

for horizontal_line=1:length(line_pos_y)
    line([-2 2], [line_pos_y line_pos_y],'color','black')
end
for pc_idx = 1:3
    xOffset = 0.25*(pc_idx-1);
    line([xOffset xOffset], [1 length(idxs)],'color','black')
end
for i=1:length(idxs)
    for pc_idx = 1:3
        pcValue = pcs(idxs(i),pc_idx);
        xOffset = 0.25*(pc_idx-1);
        if(pcValue>=0)
            rectangle('FaceColor',colors(2*(pc_idx-1)+1,:),'EdgeColor','none','Position',[xOffset,i-barWidth/2,pcValue,barWidth])
        else
            rectangle('FaceColor',colors(2*(pc_idx-1)+2,:),'EdgeColor','none','Position',[pcValue+xOffset,i-barWidth/2,-pcValue,barWidth])
        end
    end
end
ylim([1 length(idxs)]);
xlim([min(pcs(:,1))*1.1 max(pcs(:,3))*1.1+0.5]);
set(ax2,'XTick',[0, 0.25, 0.5],'XTickLabels',{'PC1', 'PC2', 'PC3'});
xlabel('Principal Components');
set(ax2,'YTick',[]);

%% subplot 3
ax3 = subplot(2,2,3);
hold on;
for pc_idx=1:3
    plot(pcs(:,pc_idx)'*DFF,'Color',colors(2*(pc_idx-1)+1,:))
end
xlim([1 size(pcs,1)])
xlabel("Time");
ylabel("Principal Component");

%% subplot 4
ax4 = subplot(2,2,4);
hold on;
bar(explained(1:10),'FaceColor',[0.5,0.5,0.5]);
plot(cumsum(explained(1:10)),'-o','color','black');
xlabel("Principal Component");
ylabel("Percent Variance Explained");
xlim([0.25 10.75]);

%% Adjust axes positions
set(ax1,'position',[.1 .5 .39 .39])
set(ax2,'position',[.5 .5 .39 .39])
set(ax3,'position',[.1 .1 .39 .39])
set(ax4,'position',[.6 .1 .3 .3])

%% PCAs over time
figure();
x=pcs(:,1)'*DFF;
y=pcs(:,2)'*DFF;
z=pcs(:,3)'*DFF;
col = (1:length(z));
surface([x;x],[y;y],[z;z],[col;col],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',2);
xlabel('PC 1');
ylabel('PC 2');
colorbar;
end
