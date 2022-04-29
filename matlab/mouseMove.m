function mouseMove(object,eventdata)
    try 
        if length(object.Children) > 1
            % for multi-plots made with subplot
            for i=1:size(object.Children)
                if( strcmp(object.Children(i).Type,'axes') )
                    C = get(object.Children(i),'CurrentPoint');
                    title(object.Children(i),['(X,Y) = (',num2str(C(1,1)),',',num2str(C(1,2)),')']);
                end
            end
        else
            % for multi-plots made with tiledlayout
            for i=1:size(object.Children.Children)
                if( strcmp(object.Children.Children(i).Type,'axes') )
                    C = get(object.Children.Children(i),'CurrentPoint');
                    title(object.Children.Children(i),['(X,Y) = (',num2str(C(1,1)),',',num2str(C(1,2)),')']);
                end
            end
        end
    catch
        % catch for single plots
        C = get(gca,'CurrentPoint');
        title(gca,['(X,Y) = (',num2str(C(1,1)),',',num2str(C(1,2)),')']);
    end
end

