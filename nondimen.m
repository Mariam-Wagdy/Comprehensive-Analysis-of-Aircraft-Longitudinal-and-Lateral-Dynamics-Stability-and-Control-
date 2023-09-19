function nondimen(y,t,c,u0)
n=400;
y(:,1)=(y(:,1)/u0)./y(:,4);
y(:,2)=(y(:,2)/u0)./y(:,4);
y(:,3)=(y(:,3)*c/(2*u0))./y(:,4);
y(:,4)=y(:,4)./y(:,4);

figure;
tiledlayout(4,1)
ax1 = nexttile;
plot(t(1:n),y(1:n,1));
title('Nondimentional x acceleration')
xlabel('time (s)')
ylabel('u/uo/theta dot')

ax2 = nexttile;
plot(t(1:n),y(1:n,2));
title('Nondimentional z acceleration')
xlabel('time (s)')
ylabel('w/uo/theta dot')

ax3 = nexttile;
plot(t(1:n),y(1:n,3));
title('Nondimentional y angular acceleration')
xlabel('time (s)')
ylabel('qc/2uo/theta dot')

ax4 = nexttile;
plot(t(1:n),y(1:n,4));
title('Nondimentional pitch angle')
xlabel('time (s)')
ylabel('theta dot/theta dot')
end
