function h = plot_MF(x, alpha, muAx)
    nInput = size(alpha,3);
    nMF = size(alpha,2);
    A = zeros(nMF, length(x), nInput);
    
    for i = 1:nInput
        for j = 1:nMF
            A(j,:,i) = muAx(x,alpha(1,j,i), alpha(2,j,i), alpha(3,j,i));
        end
    end
    
    for i = 1:nInput
        subplot(1,nInput,i);
        plot(x,A(:,:,i));
        x_txt = alpha(3,:,i);
        y_txt = 1.15*ones(size(x_txt));
        str = cell(1, nMF);
        for j = 1:nMF
            str{j} = strcat('A_{',num2str(i),num2str(j),'}');
        end
        text(x_txt, y_txt, str,'HorizontalAlignment', 'center', 'color', 'black');
        title(strcat('Membership function of x_{',num2str(i),'}'));
        axis([min(x) max(x) -0.1 1.3]);
    end
end