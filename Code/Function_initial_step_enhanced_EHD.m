function y=Function_initial_step_enhanced_EHD(P,xx,yy,a)  
y_initial=twoD(P,xx,yy);
y_second=Function_EHD2019_initial_step(P,xx,yy,a);
y=[y_initial;y_second]
y=unique(y,'rows'); 
end

