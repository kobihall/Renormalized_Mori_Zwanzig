h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
for n = 1:201
    % Draw plot for y = x.^n
    %hold on
    plot(real(u_real(:,n)))
    hold on
    plot(imag(u_real(:,n)))
    hold off
    ylim([-3 3])
    xlim([0 65])
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
      end 
  end