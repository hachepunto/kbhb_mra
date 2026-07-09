library(matrixStats)

TNBCclassif = function(x, version=c("lehmann", "bareche"), shortName=FALSE,
                       coef=FALSE, sig=NULL, rescale=TRUE)
{ version = match.arg(version);

if (rescale) { y = (x-rowMeans(x, na.rm=TRUE))/rowSds(x, na.rm=TRUE) } else { y = x; }
if (is.null(sig)) { sig = getSig("lehmann"); }
if (version=="bareche") { sig$basal_like_2 = NULL; }
if (shortName) { names(sig) = c(basal_like_1=ifelse(version=="lehmann", "BL1", "BL"), basal_like_2="BL2",
                                immunomodulatory="IM",luminal_ar="LAR",
                                mesenchymal="M",mesenchymal_stem_like="MSL")[names(sig)] }
leh = calcSig(y, sig, balanced=TRUE)
if (coef) { return(leh); }
cl = colnames(leh)[apply(leh,1,function(i) { r=which.max(i); if (length(r)==1) {r} else { NA; } }) ];
names(cl) = colnames(y);
return(cl);
}


calcSig = function(d, sig, dropEmpty=TRUE, balanced=FALSE, isCount=FALSE, useNB=FALSE,
                   updateWeights=FALSE, returnW=FALSE, maxN=200)
{ if (is(d, "Matrix")) { colMeans=Matrix::colMeans; colSums=Matrix::colSums; }
  if (is.list(sig) && !is.data.frame(sig))
  { ret = sapply(sig, function(i) calcSig(d, i, balanced=balanced, isCount=isCount, useNB=useNB, updateWeights=updateWeights));
  if (dropEmpty)
  { ret = ret[,!colSums(!is.na(ret))==0];
  }
  return(ret);
  }
  
  s = as.character(sig[,1]);
  s[s==""] = NA;
  x = intersect(rownames(d), s);
  x = x[!is.na(x)];
  if (length(x) == 0)
  { if (colnames(sig)[1] == "name") { other = "entrez"; } else { other = "name"; }
    s = as.character(sig[,other]);
    s[s==""] = NA;
    x = intersect(rownames(d), s);
    x = x[!is.na(x)];
    if (length(x) == 0)
    { warning("No gene in common");
      return(rep(NA, ncol(d)));
    }
  }
  
  sig = sig[match(x, s),];
  
  if (isCount)
  { if (useNB)
  { if (length(x)==1) { return(log10(1e4*d[x,]/colSums(d)+1)); }
    w = rowSums(d[x,])>1; x = x[w]; sig = sig[w,];
    if(!require(MASS)) { stop("MASS library not installed"); }
    if (any(duplicated(colnames(d)))) { stop("No duplicated colnames"); }
    y = cbind(as.vector(d[x,]), rep(colSums(d), each=length(x)));
    id = factor(rep(colnames(d), each=length(x)))
    gid = factor(rep(x, ncol(d)));
    a = sig$coefficient;
    Ns = colnames(d);
    if (ncol(d)>maxN)
    { Ns = colnames(d)[sample(1:ncol(d), ifelse(length(x)>maxN, maxN/2, maxN))];
    su = id %in% Ns;
    y = y[su,]; id = id[su]; gid = gid[su];
    }
    
    if (all(sig$coefficient==1))
    { fm = glm.nb.noErr(y[,1] ~ id + gid -1 + offset(log(y[,2])));
    co = coef(fm);
    co2 = co[grep("^id", names(co))];
    names(co2) = sub("^id", "", names(co2))
    co2 = co2[Ns]
    #co2 = co2[colnames(d)]; co2[is.na(co2)] = 0; names(co2) = colnames(d);
    if (updateWeights)
    { a = tapply(1:nrow(y), gid, function(i) { z=y[i,];
    coef(glm.nb.noErr(z[,1] ~ co2 + offset(log(z[,2]))))["co2"]})
    }
    } else { updateWeights=updateWeights+1; }
    
    for (iter in 1:updateWeights)
    { g2b = rep(a, length(id)/length(a))
    fm = glm.nb.noErr(y[,1] ~ g2b*id-g2b-id-1 + gid + offset(log(y[,2])));
    co = coef(fm);
    co2 = co[grep("^g2b:id", names(co))]; co2[is.na(co2)] = 0;
    names(co2) = sub("^(g2b)*:id", "", names(co2));
    co2 = co2[Ns];
    if (iter<updateWeights)
    {  a = tapply(1:nrow(y), gid, function(i) { z=y[i,];
    coef(glm.nb.noErr(glm.nb(z[,1] ~ co2 + offset(log(z[,2])))))["co2"]} );
    browser();
    }
    }
    
    if (ncol(d)>maxN)
    { co = coef(fm);
    cog = co[grep("^(gid|g2b:id)",names(co))]
    names(cog) = sub("^(gid|g2b:id)", "", names(cog))
    cog = cog[x]; cog[is.na(cog)] = 0; names(cog) = x;
    #cog = cog+co["(Intercept)"];
    if (all(a==1))
    { co2 = sapply(1:ncol(d), function(i)
    { if (sum(d[x,i])==0) { return(-Inf); }
      coef(suppressWarnings(glm(d[x,i] ~ offset(log(N[i]) + cog), family=negative.binomial(fm$theta))))
    } );
    } else
    { co2 = sapply(1:ncol(d), function(i)
    { if (sum(d[x,i])==0) { return(-Inf); }
      coef(suppressWarnings(glm(d[x,i] ~ a + offset(log(N[i]) + cog), family=negative.binomial(fm$theta))))["a"]
    } );
    }
    names(co2) = colnames(d)
    }
    
    if (returnW) { names(a) = x; return(list(sig=co2, weights=a)); }
    return(co2)
  }
    fCalc = function(d, x, coef)
    { r = colSums(d[x,,drop=FALSE]*coef, na.rm=TRUE)/colSums(d, na.rm=TRUE);
    sign(r)*log10(1+abs(r)*100);
    }
  } else { fCalc = function(d, x, coef) { colMeans(d[x,,drop=FALSE]*coef); } }
  
  if (balanced)
  { w = sig[,"coefficient"] > 0;
  #val1 = d[x,w,drop=FALSE] * sig[w,"coefficient"];
  #val1 = colMeans(val, na.rm=TRUE);
  val1 = fCalc(d, x[w], sig[w,"coefficient"]);
  #val2 = d[x,!w,drop=FALSE] * sig[!w,"coefficient"];
  #val2 = colMeans(val, na.rm=TRUE);
  val2 = fCalc(d, x[!w], sig[!w,"coefficient"]);
  val = (val1+val2)/2;
  }
  else
  { #val = d[x,,drop=FALSE] * sig[,"coefficient"];
    #val = colMeans(val, na.rm=TRUE);
    val = fCalc(d, x, sig[,"coefficient"]);
  }
  return(val);
}
                 
getSig = function(name="sig")
{ nm = paste0(baseDir(), "/Sigs/", name, ".RData");
if (!file.exists(nm)) { stop("Sig ", name, " does not exist."); }
f = load(nm);
if (!identical(f, "sig")) { stop("Not the right info in sig.RData"); }
return(sig);
}
                 
