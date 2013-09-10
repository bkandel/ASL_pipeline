#!/usr/bin/env Rscript
getPkg <- function(pckg) install.packages(pckg, repos = "http://cran.r-project.org")
pkg <- try(require(getopt))
if(!pkg) {
  cat("Installing 'getopt' from CRAN")
  getPkg('getopt')
  require('getopt')
}
spec <- matrix(c(
  'labelImage',      'l', 1, 'character', 'Name of label image', 
  'pcasl',           'p', 1, 'character', 'Name of pcasl image',
  'anatomicalImage', 'a', 1, 'character', 'Name of anatomical image',
  'gmImage',         'g', 1, 'character', 'Name of gray matter image', 
  'wmImage',         'w', 1, 'character', 'Name of white matter image', 
  'csfImage',        'c', 1, 'character', 'Name of CSF image', 
  'mask',            'm', 1, 'character', 'Name of brain mask image', 
  'dim',             'd', 1, 'integer',   'Dimensionality of images',
  'template',        't', 1, 'character', 'Name of template image', 
  'outdir',          'o', 1, 'character', 'Name of output directory', 
  'prefix',          'r', 1, 'character', 'Output prefix', 
  'tmplt2t1deform',  'e', 1, 'character', 'Template to T1 transform (deformable)', 
  'tmplt2t1affine',  'f', 1, 'character', 'Template to T1 transform (affine)'), 
               ncol=5, byrow=TRUE)
opt <- getopt(spec)
for(img in c(opt$labelImage, opt$perfusionImage, opt$anatomicalImage, 
             opt$gmImage, opt$wmImage, opt$csfImage, opt$mask, 
             opt$tmpl2t1affine, opt$tmplt2t1deform)){
  if(!file.exists(img))
    stop(paste(img, 'does not exist.'))
}
if(!file.exists(opt$outdir)){
  dir.create(opt$outdir)
}

require(ANTsR)
t1       <- antsImageRead(opt$anatomicalImage, opt$dim)
template <- antsImageRead(opt$template, opt$dim)
lab      <- antsImageRead(opt$labelImage, opt$dim)
pcasl    <- antsImageRead(opt$pcasl, (opt$dim+1))
gm       <- antsImageRead(opt$gmImage, opt$dim)
wm       <- antsImageRead(opt$wmImage, opt$dim)
csf      <- antsImageRead(opt$csfImage, opt$dim)
mask     <- antsImageRead(opt$mask, opt$dim)
outdir   <- opt$outdir

####### ASL processing.  #########
avgpcasl<-new( "antsImage" , "float" , 3 )
antsMotionCorr( list( d = 3 , a = pcasl , o = avgpcasl ) )
transform.asl2t1 <- suppressMessages(antsRegistration(t1, 
                                                      avgpcasl, "Rigid", outdir))
## should probably do something more like what we do in antsIntermodalityIntersubject.sh
concatenatedtx<-c(opt$tmplt2t1deform, opt$tmplt2t1affine,
                  unlist( transform.asl2t1$fwdtransforms ) )
print(paste('Concatenated transform list is:', concatenatedtx))
pcasl.warped <- antsApplyTransforms(fixed=template, moving=avgpcasl,
                       transformlist=concatenatedtx )
par( mfrow=c(2,1))
plotANTsImage( template, slices="60x120x10", axis=0, 
               outname=paste(outdir, '/template.png', sep='') )
plotANTsImage( pcasl.warped  , slices="60x120x10", axis=0, 
               outname=paste(outdir, '/pcasl_warped.png', sep=''))
# also get the t1 brain extraction mask into the space of the asl
pcaslmask<-new( "antsImage" , "float" , 3 )
pcaslmask <- antsApplyTransforms(fixed=avgpcasl, moving=mask,
         transformlist=transform.asl2t1$invtransforms[1], interpolator="NearestNeighbor")
par( mfrow=c(2,1))
plotANTsImage( pcaslmask, slices="1x14x2",axis=0, 
               outname=paste(outdir, '/pcasl_mask.png', sep=''))

pcasl.processing <- aslPerfusion(pcasl, mask=pcaslmask, moreaccurate=FALSE )
# in above, we should truncate outliers or add a 2nd filter / function to do so
pcasl.perfusion <- pcasl.processing$perfusion

pcasl.parameters <- list( sequence="pcasl", m0=pcasl.processing$m0 )
cbf <- quantifyCBF( pcasl.perfusion, pcaslmask, pcasl.parameters )
cbf.warped2template <- antsApplyTransforms( template, cbf$meancbf,
                         transformlist=concatenatedtx )
cbf.warped2t1 <- antsApplyTransforms( fixed=t1, 
                moving=cbf$meancbf, 
                transformlist=c(unlist(transform.asl2t1$fwdtransforms ) ))

gm.mask <- new('antsImage', 'float', 3)
wm.mask <- new('antsImage', 'float', 3)
ThresholdImage(3, gm, gm.mask, 0.5, 999)
ThresholdImage(3, wm, wm.mask, 0.5, 999)
cbfvals.gm <- getROIValues(cbf.warped2template, lab, gm.mask)
cbfvals.wm <- getROIValues(cbf.warped2template, lab, wm.mask)

####### Write out images. #####
antsImageWrite(cbf.warped2template, 
               paste(outdir, '/', opt$prefix, 'MeanCBFWarpedToTemplate.nii.gz',
                                          sep=''))
antsImageWrite(cbf.warped2t1, 
    paste(outdir, '/', opt$prefix, 'MeanCBFWarpedToT1.nii.gz', sep=''))
antsImageWrite(cbf$meancbf, paste(outdir, '/', opt$prefix, 'MeanCBF.nii.gz', sep=''))
antsImageWrite(avgpcasl, paste(outdir, '/', opt$prefix, 'AveragePCASL.nii.gz', sep=''))
myvals <- data.frame(GMVals=cbfvals.gm$roiMeans, WMVals=cbfvals.wm$roiMeans)
row.names(myvals) <- as.character(cbfvals.gm$roiValues)
write.csv(myvals, paste(outdir, '/', opt$prefix, 'ROIValues.csv', sep=''))
