clear all
clc

results = [];

sections = [];

data = importdata('output_Epfield_xinc.dat');

tamanho=length(data)

dataneg=[];
datapos=[];
NormEback=[];
NormEforw=[];



for i =1:(tamanho/2)
    dataneg(i,:)=data(2*i-1,:);
    datapos(i,:)=data(2*i,:);
end

Eback=[];
Eforw=[];
Eback=[dataneg(:,5)+1j*dataneg(:,6),dataneg(:,7)+1j*dataneg(:,8),dataneg(:,9)+1j*dataneg(:,10)];
Eforw=[datapos(:,5)+1j*datapos(:,6),datapos(:,7)+1j*datapos(:,8),datapos(:,9)+1j*datapos(:,10)];

for i = 1:tamanho/2
  NormEback(i)=norm(Eback(i,:)); 
  NormEforw(i)=norm(Eforw(i,:)); 
end
plot(dataneg(:,1), NormEback)
%plot(datapos(:,1), NormEforw)

%plot(datapos(:,1), sqrt(datapos(:,5).^2+datapos(:,6).^2))
