x_plot = linspace(0,1,1000)
a_y = 100
b_y = 1000
y_plot = gampdf(x_plot,a_y,1/b_y)
figure(1)
plot(x_plot, y_plot)