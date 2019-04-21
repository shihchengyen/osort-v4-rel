function changeBARcolor(childrenhandles,barcolor)


for i=1:length(childrenhandles)
    try
       set(childrenhandles,'FaceColor',barcolor);
        set(childrenhandles,'EdgeColor',[0 0 0]);
    end
end