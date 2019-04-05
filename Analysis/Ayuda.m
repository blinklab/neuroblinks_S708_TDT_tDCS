idx= find(trials.t_ONOFF & trials.t_polarity) %busca en ambas columnas, los valores que sean mayores de 1
idx2= find(trials.t_ONOFF & ~trials.t_polarity) %busca en la primera columna los valores que sean mayores de 1 y en la 2 los que tengan valor 0
idx= find(trials.t_ONOFF & trials.t_polarity,5,'last') %me devuelve solo los 5 ultimos
idx= find(trials.t_ONOFF & trials.t_polarity & trials.session_of_day==3) %para coger solo los de la sesion 3
idxT= [idx, idx2] %creamos una matriz con esas dos columnas
idxC= find(~trials.t_ONOFF & ~trials.t_polarity,5,'last') %seleccionamos control
idxr=find(trials.session_of_day==3 & ismember(trials.trialnum, idx)) %Me da las coordenadas que cumplen que es la sesion3 y que el numero de trial son los numeros de la matriz idx (que son los trialnum donde aplico anodal, por ejemplo)


plot(trials.tm(1,:),trials.eyelidpos(idxC,:),'k') %figura control color negro
figure
plot(trials.tm(1,:),trials.eyelidpos(idx,:)) %figura catodal
figure
plot(trials.tm(1,:),trials.eyelidpos(idx2,:)) %figura anodal
hold on %mantenemos la primera figura y le agregamos la siguiente encima

figure
plot(trials.tm(1,:),trials.eyelidpos(intersect(trials.trialnum,idx),:))
title('catodal')
ylabel('eyelidpos')
xlabel('time')
prueba(3).Color = 'r' %al trazo 3 de la figura 'prueba' le ponemos color rojo


%Codigo para representar figura con trialddata abierto
idx= find(trials.t_ONOFF & trials.t_polarity & trials.session_of_day==3) %para coger solo los de la sesion 3
idx2= find(trials.t_ONOFF & ~trials.t_polarity & trials.session_of_day==3) %para coger solo los de la sesion 3
idxC= find(~trials.t_ONOFF & ~trials.t_polarity,5,'last') %seleccionamos control
plot(trials.tm(1,:),trials.eyelidpos(idxC,:),'k') %figura control color negro
hold on
plot(trials.tm(1,:),trials.eyelidpos(idx,:),'B') %figura control color negro
plot(trials.tm(1,:),trials.eyelidpos(idx2,:),'r') %figura control color negro
