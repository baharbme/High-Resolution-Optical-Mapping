function checkresult_d(t_data, x, dMin, tdMin, dMax, tdMax)
scrsz = get(0,'ScreenSize');
figure('Position',[100 scrsz(4)/3 scrsz(3)/2 scrsz(4)/2])
scatter(t_data,x,'linewidth',2);
hold on; plot([tdMin, tdMax],[dMin,dMax],'g*')
box on; grid on
xlabel('time (ms)'); ylabel('intensity')