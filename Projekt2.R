library("TSdist")
library("ggplot2")
library("dplyr")
library("tidyr")
library("psych")
library("stats")
library("cluster")
library("kableExtra")
library("gridExtra")

#========
piwo<-read.csv("E:/Studia/Sem. 5/SADzR/Projekt2/piwo.csv",sep = ";", dec = ".")
rownames(piwo)<-piwo$marka
piwo<-piwo[,-1]

#View(piwo)
#dist(piwo)
#mahalanobis(as.numeric(piwo[1,]),as.numeric(piwo[2,]),var(piwo))
#EuclideanDistance(as.numeric(piwo[1,]),as.numeric(piwo[2,]))
#TSDistances(as.numeric(piwo[1,]),as.numeric(piwo[2,]),distance ="dissim")
#EuclideanDistance(c(1,2,3,4,1,0),c(1,2,3,4,0,1))
piwo<-rbind(piwo, wzor=apply(piwo_transform,2,max))
dist_vect<-as.vector(daisy(piwo,metric ="gower"))
dist_mat<-as.matrix(dist(piwo))
dist_vect2<-dist_mat[nrow(dist_mat),1:(nrow(dist_mat)-1)]
dist_vect2[2]




View(dist(piwo))
dist_mat<-as.matrix(dist_vect)

piwo_transform<-piwo
piwo_transform$zawartosc.alk<--abs(6-piwo_transform$zawartosc.alk)
piwo_transform$cena<--piwo_transform$cena

piwo_transform<-as.data.frame(scale(piwo_transform))
piwo_transform<-piwo_transform*as.list(c(5,1,1,1,1))




piwo_transform_max<-apply(piwo_transform,2,max)
wynik_vect<-NULL
for (i in 1:nrow(piwo_transform)) {
  wynik<-EuclideanDistance(as.numeric(piwo_transform[i,]),piwo_transform_max)
  wynik_vect<-c(wynik_vect,wynik)
}
piwo_transform<-cbind(piwo_transform,EDist=wynik_vect)

#View(piwo_transform)

helwig_piwo<-1-(piwo_transform[,ncol(piwo_transform)]/(mean(piwo_transform[,ncol(piwo_transform)]) + 2 * sd(piwo_transform[,ncol(piwo_transform)])))
piwo_transform<-cbind(piwo_transform,helwig = helwig_piwo)

piwo_rank<-cbind(piwo,helwig = helwig_piwo)
piwo_rank<-cbind(piwo_rank[order(-piwo_rank$helwig),],rank=1:nrow(piwo_transform))

View(piwo_rank)
#========

# funkcja ranking nadaje ranguje obserwacje za pomoc¹ 3 dostêpnych metod
# dframe - ramka danych
# method - metoda("helwig","stand_sum", "suma_rang"), domyœlnie helwig
# var_types_vect - wektor zawieraj¹cy charakter zmiennych znajduj¹cych siê 
#   w ka¿dej kolumnie(s-stymulanta,d-destymulanta,n-nominanta), domyœlnie wszystkie stymulanty
# optimum_vect - wektor wartoœci optymalnych
# weight_vect - wektor wag dla zmiennych w ka¿dej kolumnie,domyœlnie wszystkie wagi jednakowe
# weight_auto - wartoœæ logiczna(TRUE/FALSE), okreœlaj¹ca czy wagi zmiennych maj¹ zostaæ wygenerowane
#   automatycznie, domyœlnie FALSE

ranking<-function(dframe, method="helwig", var_types_vect = rep("s",ncol(dframe)), optimum_vect, 
                  weight_vect = rep(1,ncol(dframe)),weight_auto=FALSE){
  
  dframe_transform<-dframe
  i_opt<-0
  #zamiana nominant i destymulant na stymulanty
  for (i_col in 1:ncol(dframe)) {
    if (var_types_vect[i_col]=="d") {
      dframe_transform[,i_col]<--dframe_transform[,i_col]
    } 
    if (var_types_vect[i_col]=="n") {
      i_opt<-i_opt+1
      dframe_transform[,i_col]<--abs(optimum_vect[i_opt]-dframe_transform[,i_col])
    }
  }
  #skalowanie wartosci zmiennych
  dframe_transform<-as.data.frame(scale(dframe_transform))
  
  #nadawanie wag zmiennych, w przypadku wyboru opcji automatycznego doboru
  if (weight_auto==TRUE) {
    sd_mean<-apply(dframe_transform,2,sd)/apply(dframe_transform,2,mean)
    weight_vect<-sd_mean/sum(sd_mean)
  }
  
  #porzadkowanie metoda helwiga
  if(method=="helwig"){
    #nalozenie wag na standaryzowane dane
    dframe_transform<-dframe_transform*as.list(weight_vect)
    #ustalenie wektora zawieraj¹cego najlepsze wartosci("wzor")
    dframe_transform_max<-apply(dframe_transform,2,max)
    distance_vect<-NULL
    #obliczenie odleglosci euklidesowych kazdej z oberwacji od "wzoru"
    for (i in 1:nrow(dframe_transform)) {
      distance<-EuclideanDistance(as.numeric(dframe_transform[i,]),dframe_transform_max)
      distance_vect<-c(distance_vect,distance)
    }
    dframe_transform<-cbind(dframe_transform,EDist=distance_vect)
    #wyznaczenie miary dla kazdego z obiektow
    helwig_dframe<-1-(dframe_transform[,ncol(dframe_transform)]/(mean(dframe_transform[,ncol(dframe_transform)]) + 2 * sd(dframe_transform[,ncol(dframe_transform)])))
    dframe_transform<-cbind(dframe_transform,helwig = -helwig_dframe)
    
    #dopisanie pozycji do kazdego obiektu
    dframe_unordered<-cbind(dframe,helwig = helwig_dframe,rank=round(rank(-helwig_dframe, ties.method = "first"),digits = 0))
    #uprzadkowanie wg. rankingu
    dframe_rank<-dframe_unordered[order(dframe_unordered$rank),]
    
  }
  #porzadkowanie metoda standaryzowanych sum
  if(method=="stand_sum"){
    #nalozenie wag na standaryzowane dane
    dframe_transform<-dframe_transform*as.list(weight_vect)
    
    #uzyskanie wektora sum dla kazdej obserwacji, a nastepnie standaryzacja sum
    stand_sum_vect<-apply(dframe_transform,1,sum)
    stand_sum_vect<-(stand_sum_vect-min(stand_sum_vect))/(max(stand_sum_vect-min(stand_sum_vect)))
    dframe_transform<-cbind(dframe_transform,stand_sum = stand_sum_vect)
    
    #dopisanie pozycji do kazdego obiektu
    dframe_unordered<-cbind(dframe,stand_sum = stand_sum_vect,rank=round(rank(-stand_sum_vect, ties.method = "first"),digits = 0))
    #uprzadkowanie wg. rankingu
    dframe_rank<-dframe_unordered[order(dframe_unordered$rank),]
  }
  #porzadkowanie metoda sumy rang
  if(method=="suma_rang"){
    #utworzenie dodatkowych kolumn zawierajacych rangi wg. kazdej ze zmiennych
    for (i_col in 1:ncol(dframe)) {
      new_col<-paste0("rank",i_col)
      dframe_transform<-cbind(dframe_transform,rank(-dframe_transform[,i_col],ties.method = "average"))
      colnames(dframe_transform)[i_col+ncol(dframe)]<-new_col
    }
    #nalozenie wag na rangi
    dframe_transform[,(ncol(dframe)+1):ncol(dframe_transform)]<-dframe_transform[,(ncol(dframe)+1):ncol(dframe_transform)]*as.list(weight_vect)
    
    #zsumowanie rang 
    rank_sum_vect<-apply(dframe_transform[,(ncol(dframe)+1):ncol(dframe_transform)],1,mean)
    dframe_transform<-cbind(dframe_transform,rank_sum = -rank_sum_vect)
    
    #dopisanie pozycji do kazdego obiektu
    dframe_unordered<-cbind(dframe,rank_sum = rank_sum_vect,rank=round(rank(rank_sum_vect, ties.method = "first"),digits = 0))
    #uprzadkowanie wg. rankingu
    dframe_rank<-dframe_unordered[order(dframe_unordered$rank),]
  }
  #zwrocenie wejsciowej ramki danych z dodanymi rangami dla kazdej obserwacji
  return(list(dframe_unordered,dframe_rank,dframe_transform))  
}

#=====================================
sr_df<-(ranking(piwo, "suma_rang", c("n","d","s","s","s"), c(6),c(1,1,1,1,1)))
ss_df<-(ranking(piwo, "stand_sum", c("n","d","s","s","s"), c(6),c(1,1,1,1,1)))
h_df<-ranking(piwo, "helwig", c("n","d","s","s","s"), c(6),c(1,1,1,1,1))
rank_df<-data.frame(nazwa=rownames(sr_df),suma_rang=as.numeric(sr_df$rank),stand_sum=as.numeric(ss_df$rank),helwig=as.numeric(h_df$rank))

ggplot(data=rank_df,aes(x=suma_rang,y=nazwa),alpha=0.02)+
  geom_point(aes(x=suma_rang,y=nazwa),col="blue",size=4,alpha=0.4,position=position_jitter(height = 0.4,width=0))+
  geom_point(aes(x=stand_sum,y=nazwa),col="red",size=4,alpha=0.4,position=position_jitter(height = 0.4,width=0))+
  geom_point(aes(x=helwig,y=nazwa),col="green",size=4,alpha=0.4,position=position_jitter(height = 0.4,width=0))+
  scale_x_discrete(limits=1:20)



#dane-kultura
dane<-read.csv("E:/Studia/Sem. 5/SADzR/Projekt2/Dane_kultura/dane_kultura.csv",sep = ";", dec = ",")
dane<-na.omit(dane)

dane_wsk<-dane
dane_wsk[,3:11]<-dane[,3:11]*10000/dane$ludnosc
dane_wsk<-dane_wsk[rank(-dane_wsk$ludnosc)<=25,]

rownames(dane_wsk)<-dane_wsk$Nazwa
dane_wsk<-dane_wsk[,-c(1:2)]

rownames(dane)<-dane$Nazwa
dane<-dane[,-c(1)]

View(dane)
View(dane_wsk)

View(ranking(dane_wsk, "stand_sum",c("s","s","s","s","s","s","s","s","d")))

#=====================================

#dane_ekologia
dane<-read.csv("E:/Studia/Sem. 5/SADzR/Projekt2/Dane_ekologia/dane_ekologia.csv",sep = ";", dec = ",")
dane<-dane[,1:10]
dane<-na.omit(dane)
View(dane)

dane_wsk<-dane
dane_wsk[,c(4:6,10)]<-dane[,c(4:6,10)]*1000/dane$ludnosc
dane_wsk<-dane_wsk[rank(-dane_wsk$ludnosc)<=25,]

rownames(dane_wsk)<-dane_wsk$Nazwa
dane_wsk<-dane_wsk[,-c(1:2)]
View(dane_wsk)

r20<-(ranking(dane_wsk, "helwig",c("s","d","d","d","s","n","d","s"),c(6)))
r25<-(ranking(dane_wsk, "helwig",c("s","d","d","d","s","n","d","s"),c(6)))

View(r20)
View(r25)

kable((ranking(dane_wsk, "helwig",c("s","d","d","d","s","n","d","s"),c(6)))) %>% kable_styling(full_width = F)






w<-c()
for (i in 1:ncol(dane_wsk)) {
  breaks <- pretty(range(dane_wsk[,i]), n = nclass.FD(dane_wsk[,i]), min.n = 1)
  bwidth <- breaks[2]-breaks[1]
  w[i]<-bwidth
}

w1<-dane_wsk %>%
  ggplot(aes(udzial_terenow_ziel)) +
  geom_histogram(binwidth=w[1],fill="red")+
  labs(y="")
w2<-dane_wsk %>%
  ggplot(aes(pojazdy_paliwa_plynne)) +
  geom_histogram(binwidth=w[2],fill="orange")+
  labs(y="")
w3<-dane_wsk %>%
  ggplot(aes(zuzycie_wody)) +
  geom_histogram(binwidth=w[3],fill="yellow")+
  labs(y="")
w4<-dane_wsk %>%
  ggplot(aes(odpady_niezneutralizowane)) +
  geom_histogram(binwidth=w[4],fill="green")+
  labs(y="")
w5<-dane_wsk %>%
  ggplot(aes(odpady_selektywnie_do_ogolu)) +
  geom_histogram(binwidth=w[5],fill="blue")+
  labs(y="")
w6<-dane_wsk %>%
  ggplot(aes(wsk_ladunek_scieki)) +
  geom_histogram(binwidth=w[6],fill="purple")+
  labs(y="")
w7<-dane_wsk %>%
  ggplot(aes(wsk_jakosc_powietrza)) +
  geom_histogram(binwidth=w[7],fill="pink")+
  labs(y="")
w8<-dane_wsk %>%
  ggplot(aes(wydatki_ochrona_srod)) +
  geom_histogram(binwidth=w[8],fill="white")+
  labs(y="")

grid.arrange(w1, w2, w3, w4,w5,w6,w7,w8, 
             ncol = 2, nrow = 4,top="Histogramy dla zmiennych charakteryzuj¹cych miasta") 




w1<-cbind(dane_wsk,czynnik=rep(colnames(dane_wsk)[1],nrow(dane_wsk)))    %>%
  ggplot(aes(x = czynnik, y = udzial_terenow_ziel)) +
  geom_boxplot(fill="red")+
  labs(x="")
w2<-cbind(dane_wsk,czynnik=rep(colnames(dane_wsk)[2],nrow(dane_wsk)))    %>%
  ggplot(aes(x = czynnik, y = pojazdy_paliwa_plynne)) +
  geom_boxplot(fill="orange")+
  labs(x="") 
w3<-cbind(dane_wsk,czynnik=rep(colnames(dane_wsk)[3],nrow(dane_wsk)))    %>%
  ggplot(aes(x = czynnik, y = zuzycie_wody)) +
  geom_boxplot(fill="yellow")+
  labs(x="") 
w4<-cbind(dane_wsk,czynnik=rep(colnames(dane_wsk)[4],nrow(dane_wsk)))    %>%
  ggplot(aes(x = czynnik, y = odpady_niezneutralizowane)) +
  geom_boxplot(fill="green")+
  labs(x="")
w5<-cbind(dane_wsk,czynnik=rep(colnames(dane_wsk)[5],nrow(dane_wsk)))    %>%
  ggplot(aes(x = czynnik, y = odpady_selektywnie_do_ogolu)) +
  geom_boxplot(fill="blue")+
  labs(x="") 
w6<-cbind(dane_wsk,czynnik=rep(colnames(dane_wsk)[6],nrow(dane_wsk)))    %>%
  ggplot(aes(x = czynnik, y = wsk_ladunek_scieki)) +
  geom_boxplot(fill="purple")+
  labs(x="") 
w7<-cbind(dane_wsk,czynnik=rep(colnames(dane_wsk)[7],nrow(dane_wsk)))    %>%
  ggplot(aes(x = czynnik, y = wsk_jakosc_powietrza)) +
  geom_boxplot(fill="pink")+
  labs(x="") +
  theme_bw()
w8<-cbind(dane_wsk,czynnik=rep(colnames(dane_wsk)[8],nrow(dane_wsk)))    %>%
  ggplot(aes(x = czynnik, y = wydatki_ochrona_srod)) +
  geom_boxplot(fill="black")+
  labs(x="")+
  theme_bw()

grid.arrange(w1, w2, w3, w4,w5,w6,w7,w8, 
             ncol = 4, nrow = 2,top="Wykresy pude³kowe dla zmiennych charakteryzuj¹cych miasta")                           

kable(describe(dane_wsk)[c(2:5,8:13)],digits=2, align = "c") %>% kable_styling(full_width = F,position="center",bootstrap_options = c("striped", "hover", "condensed", "responsive"))
kable(dane_wsk,digits=2, align = "c") %>% kable_styling(full_width = F,position="center",bootstrap_options = c("striped", "hover", "condensed", "responsive"))


kable(dane_wsk,digits=2, align = "c") %>% kable_styling(full_width = F,position="center",bootstrap_options = c("striped", "hover", "condensed", "responsive"))%>%
  column_spec(c(2,6,9), bold = T, color = "black", background = "lightgreen")%>%
  column_spec(c(3:5,8), bold = T, color = "black", background = "#F08080")%>%
  column_spec(c(7), bold = T, color = "black", background = "lightblue")


round((ranking(dane_wsk, "helwig",c("s","d","d","d","s","n","d","s"),c(6)))[[3]],digits=2) %>%
  mutate_if(is.numeric, function(x) {
    cell_spec(x, bold = T, 
              color = spec_color(x, end = 1),
              font_size = spec_font_size(x))
  }) %>%
  kable(escape=F,digits=0, align = "c") %>% kable_styling(full_width = F,position="center",bootstrap_options = c("striped", "hover", "condensed", "responsive"))


library(maps)
library(SmarterPoland)

#dane_ekologia
dane<-read.csv("E:/Studia/Sem. 5/SADzR/Projekt2/Dane_ekologia/dane_ekologia.csv",sep = ";", dec = ",")
dane<-dane[,1:10]
dane<-na.omit(dane)
View(dane)

dane_wsk<-dane
dane_wsk[,c(4:6,10)]<-dane[,c(4:6,10)]*1000/dane$ludnosc
dane_wsk<-dane_wsk[rank(-dane_wsk$ludnosc)<=25,]

rownames(dane_wsk)<-dane_wsk$Nazwa
dane_wsk<-dane_wsk[,-c(1:2)]


poland <- map("world","poland", fill=T, col="white")
miasta<-getMPpowiaty()
data(world.cities)
cities_lon_lat <- world.cities[!duplicated(world.cities$name),]
rownames(cities_lon_lat) = cities_lon_lat[,1]
cities_lon_lat <- cities_lon_lat[cities_lon_lat$country.etc == "Poland",]
nazwy_miast<-rownames(dane_wsk)
nazwy_miast[7]<-"Cracow"
nazwy_miast[9]<-"Warsaw"
nazwy_miast[13]<-"Bielsko-Biala"
lon_miasta<-(sapply(nazwy_miast,function(x){cities_lon_lat[cities_lon_lat$name == x,]$lon}))
lat_miasta<-(sapply(nazwy_miast,function(x){cities_lon_lat[cities_lon_lat$name == x,]$lat}))
cities_coords<-data.frame("miasto"=rownames(dane_wsk),"lon"=unlist(lon_miasta),"lat"=unlist(lat_miasta),stringsAsFactors = FALSE)
rownames(cities_coords)<-NULL
poland<-data.frame("lon"=poland$x,"lat"=poland$y)

ggplot() + geom_polygon(data = poland, aes(x=lon, y = lat),fill="#ffb3b3",col="black") + 
  coord_fixed(1.3)+
  geom_point(data=cities_coords, aes(x=lon, y = lat), col="black")+
  theme_nothing()+
  geom_text(data = cities_coords[-c(14,16,17,20,22),], aes(x=lon, y = lat,label = miasto), angle=20,vjust = 0, nudge_y = -0.2, size=3,color = "black")+
  geom_text(data = cities_coords[cities_coords$miasto=="Gliwice",], aes(x=lon, y = lat,label = miasto), angle=20,hjust = 1.1, vjust = 0.2, size=3,color = "black")+
  geom_text(data = cities_coords[cities_coords$miasto=="Sosnowiec",], aes(x=lon, y = lat,label = miasto), angle=20,hjust = -0.2, vjust = 0.2, size=3,color = "black")+
  geom_text(data = cities_coords[cities_coords$miasto=="Katowice",], aes(x=lon, y = lat,label = miasto), angle=20,hjust = -0.2, nudge_y = -0.13, size=3,color = "black")+
  geom_text(data = cities_coords[cities_coords$miasto=="Bytom",], aes(x=lon, y = lat,label = miasto), angle=20,hjust = -0.2, nudge_y = -0.05, size=3,color = "black")+
  geom_text(data = cities_coords[cities_coords$miasto=="Zabrze",], aes(x=lon, y = lat,label = miasto), angle=20,hjust = 0, nudge_y = -0.08, size=3,color = "black")
  



