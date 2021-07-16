function CorrViewer(cMat)

currentPixel = [1,1];
h_fig = figure('WindowButtonDownFcn', @ChangePtPos);
h_ax = axes('Parent', h_fig);
imagesc(h_ax, reshape(cMat(1,:), 256, 256), [0 1]);
axis(h_ax, 'off', 'image');
hold(h_ax, 'on');
plot(h_ax, currentPixel(1), currentPixel(2), 'or');
hold(h_ax, 'off');

    function ChangePtPos(Obj, Evnt)
        Pos = round(Obj.Children(1).CurrentPoint);
        
        Pos = Pos(1,1:2);
        if( any(Pos < 1) | any(Pos > 256) )
            return;
        end
        currentPixel = Pos;
                
        Id = floor((currentPixel(1)-1))*256 + floor(currentPixel(2));
        if( Id < 1 )
            Id = 1;
        end
        if( Id > 256*256)
            Id = 256*256;
        end
        imagesc(h_ax, reshape(cMat(:,Id), 256, 256), [0 1]);
        axis(h_ax, 'off', 'image');
        hold(h_ax, 'on');
        plot(h_ax, currentPixel(1), currentPixel(2), 'or');
        hold(h_ax, 'off');
    end
end