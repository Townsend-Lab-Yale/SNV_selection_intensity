
# function that nudges data labels apart for plotting

data.push.function <- function(data.to.push,
                               x.min,
                               x.max,
                               touching.distance,
                               pushing.distance,
                               gamma.min = 1,
                               x.data = 'gamma_epistasis',
                               cat.data = 'tumor_type',
                               freq.data = 'Prop_tumors_with_specific_mut',
                               max.iter = 1e4){
  
  data.to.push$new_x <- data.to.push[,x.data]
  data.to.push <- data.to.push[which(data.to.push[,colnames(data.to.push)[which(colnames(data.to.push)==x.data)]]>gamma.min),]
  # data.to.push <- subset(data.to.push, x.data>gamma.min)
  
  categories <- (unique(data.to.push[cat.data]))
  
  
  for(i in 1:nrow((unique(data.to.push[cat.data])))){ #For each category of data
    still.pushing <- T;this.iter <- 1 #Breaks for the while loop 
    this.data <- data.to.push[which(data.to.push[cat.data]==as.character(categories[i,])),x.data] #just extract the data in this category
    this.order <- order(this.data)
    if(length(this.data)>1){
      #While loop to do the pushing. 
      while(still.pushing==T & this.iter <= max.iter){ 
        push.this.round <- F
        
        #loop that detects neighbors
        for(q in 1:length(this.data)){
          j <- this.order[q]
          if(length(which(abs(this.data[j]-this.data)[-j] < touching.distance))>0){ #if there is a neighbor too close
            push.this.round <- T #a push occurred 
            # need to remove the individual responsible for the pushing
            self.position <- which(which(abs(this.data[j]-this.data) < touching.distance)==j) 
            neighbors.to.push <- which(abs(this.data[j]-this.data) < touching.distance)[-self.position]
            
            #for all touching neighbors, push them both away 
            #TODO: right now data gets pushed out of the range of the data
            #      Doesn't break anything, but it would be nice to fix in case
            #      we want to constrain the spread of push. 
            for(k in 1:length(neighbors.to.push)){
              
              if(this.data[j] <= this.data[neighbors.to.push[k]]){
                if( (this.data[j] - pushing.distance) > x.min){
                  this.data[j] <- this.data[j] - pushing.distance
                  this.data[neighbors.to.push[k]] <- this.data[neighbors.to.push[k]] + pushing.distance
                }else{
                  this.data[neighbors.to.push[k]] <- this.data[neighbors.to.push[k]] + pushing.distance
                }
              }
              if(this.data[j] >= this.data[neighbors.to.push[k]]){
                if( (this.data[j] + pushing.distance) < x.max){
                  this.data[j] <- this.data[j] + pushing.distance
                  this.data[neighbors.to.push[k]] <- this.data[neighbors.to.push[k]]-pushing.distance
                }else{
                  this.data[neighbors.to.push[k]] <- this.data[neighbors.to.push[k]]-pushing.distance
                }
              }
              # if(this.data[j] == neighbors.to.push[k]){
              #   
              # }
              
            }
            
          }
          if(this.data[j] < x.min){
            push.this.round <- T #a push occurred 
            this.data[j] <- this.data[j] + pushing.distance 
          }
          if(this.data[j] > x.max){
            push.this.round <- T #a push occurred 
            this.data[j] <- this.data[j] - pushing.distance 
          }
          
        }
        
        
        this.iter <- this.iter+1
        if(this.iter==max.iter){print(paste("Maximum iterations reached! Category: ",i,sep=""))}      
        if(!push.this.round){still.pushing <- F} #if nothing happened, time to move on
      }
      data.to.push[which(data.to.push[cat.data]==as.character(categories[i,])),"new_x"] <- this.data
    }else{
      data.to.push[which(data.to.push[cat.data]==as.character(categories[i,])),"new_x"] <- this.data
    }
    message(paste("Finished category ",i," of ",nrow(categories),".",sep=""))
  }
  data.to.push$col <- NA
  for(i in 1:nrow(data.to.push)){
    if(data.to.push[i,freq.data]<0.01){
      data.to.push$col[i] <- "<1%"
    }else{
      if(data.to.push[i,freq.data]<0.02){
        data.to.push$col[i] <- "1–2%"
      }else{
        if(data.to.push[i,freq.data]<0.03){
          data.to.push$col[i] <- "2–3%"
        }else{
          if(data.to.push[i,freq.data]<0.05){
            data.to.push$col[i] <- "3–5%"
          }else{
            if(data.to.push[i,freq.data]<0.1){
              data.to.push$col[i] <- "5–10%"
            }else{
              data.to.push$col[i] <- ">10%"
            }
          }
        }
      }
    }
  }
  
  data.to.push$col <- factor(data.to.push$col, levels = rev(c("<1%","1–2%","2–3%","3–5%","5–10%",">10%")))
  return(data.to.push)
}

