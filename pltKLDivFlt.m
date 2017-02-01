clear all;close all;

comp = [0.179919282	0.274000641	0.341709108	0.321169274	0.35964163	0.22699801	0.191709837	0.189772126	0.188387062	0.187472174	0.186733606
0.179919282	0.156451867	0.100414074	0.077348173	0.068434293	0.064547694	0.063128322	0.060863267	0.059391989	0.058324609	0.057604721
0.179919282	0.099731186	0.063402162	0.040151873	0.028237542	0.021538216	0.017461222	0.014940563	0.013192255	0.011946211	0.011049
];

% 
% 
% comp = [0.179919282	0.443459623	0.292742999	0.233375098	0.228326735	0.244712108	0.254597582	0.261975656	0.267882775	0.271767932
% 0.179919282	0.181269471	0.142279152	0.120289616	0.107621729	0.100683222	0.096332814	0.093130671	0.091960885	0.088777745
% ];

% h00 = figure('position',[100 100 600 600],'Color',[1 1 1]);
% plot(0:size(comp,2)-1,comp(1,:),'s--');hold on;
% plot(0:size(comp,2)-1,comp(2,:),'r*-');hold on;
% legend('k=1','k=10');
% xlabel('time step');
% ylabel('K-L divergence');
% set(gca,'FontSize',14);
% 
% 

% 'g',x,y2,'b--o',x,y3,'c*'

h00 = figure('position',[100 100 600 600],'Color',[1 1 1]);
plot(0:size(comp,2)-1,comp(1,:),'s-');hold on;
plot(0:size(comp,2)-1,comp(2,:),'r--d');hold on;
plot(0:size(comp,2)-1,comp(3,:),'b*-');hold on;

legend('k=1','k=2','k=10');
xlabel('time step');
ylabel('K-L divergence');
set(gca,'FontSize',14);
% 
% 
% h01 = figure('position',[100 100 600 600],'Color',[1 1 1]);
% plot(0:size(comp2,2)-1,comp2(1,:),'s-');hold on;
% plot(0:size(comp2,2)-1,comp2(2,:),'g--d');hold on;
% plot(0:size(comp2,2)-1,comp2(3,:),'c*-');hold on;
% legend('k=1(fixed)','k=2(fixed)','k=10(fixed)');
% xlabel('time step');
% ylabel('K-L divergence');
% set(gca,'FontSize',14);