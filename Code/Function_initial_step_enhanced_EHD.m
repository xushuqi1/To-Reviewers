function y=Function_initial_step_enhanced_EHD(P,xx,yy,a)  
y_initial=twoD(P,xx,yy);
y_second=EHD2019(P,xx,yy,a);
y=[y_initial;y_second]
y=unique(y,'rows'); 
end

