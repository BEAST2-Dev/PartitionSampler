package beast.evolution.likelihood;

import beast.app.BeastMCMC;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.ChangeType;
import beast.evolution.alignment.GeneralUnitAlignment;
import beast.evolution.sitemodel.DPSiteModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Tree;
import sun.java2d.loops.FillRect;

/**
 * @author Chieh-Hsi Wu
 */
public class GeneralUnitDPTreeLikelihood extends DPTreeLikelihood {


    public void initAndValidate() throws Exception{
        useThreads = useThreadsInput.get() && (BeastMCMC.m_nThreads > 1);
        useThreadsEvenly = useThreadsEvenlyInput.get() && (BeastMCMC.m_nThreads > 1);
        dpVal = dpValInput.get();
        if(!(siteModelInput.get() instanceof DPSiteModel)){
            throw new RuntimeException("DPSiteModel required for site model.");
        }
        dpSiteModel = (DPSiteModel) siteModelInput.get();


        alignment = dataInput.get();
        int patternCount = alignment.getPatternCount();



        int[][] clusterWeights = new int[dpSiteModel.getSiteModelCount()][patternCount];

        int siteModelCount = dpSiteModel.getSiteModelCount();

        int siteCount = alignment.getSiteCount();
        GeneralUnitAlignment data = (GeneralUnitAlignment)alignment;
        for(int i = 0; i < siteCount; i++){
            clusterWeights[dpVal.getCurrCategory(data.getUnitBySite(i))][alignment.getPatternIndex(i)]++;
        }

        for(int i = 0; i < siteModelCount;i++){
            NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(
                    clusterWeights[i],
                    alignment,
                    (Tree) treeInput.get(),
                    useAmbiguitiesInput.get(),
                    dpSiteModel.getSiteModel(i),
                    branchRateModelInput.get());
            treeLiks.add(treeLik);

        }
        if(!(alignment instanceof GeneralUnitAlignment)){
            throw new RuntimeException("Must use GeneralUnitAlignment object.");
        }



    }

    public void splitTreeLikelihood(){

        SiteModel siteModel = dpSiteModel.getSiteModel(dpSiteModel.getLastAddedIndex());
        int[] clusterUnits = dpVal.getClusterSites(dpSiteModel.getLastAddedIndex());


        int[] patternWeights = new int[alignment.getPatternCount()];

        int prevCluster;
        if(changeType == ChangeType.SPLIT || changeType == ChangeType.SPLIT_AND_VALUE_CHANGE){
            prevCluster = dpSiteModel.getDirtySiteModelIndex();
        }else if(changeType == ChangeType.ADDED){
            prevCluster = dpSiteModel.getPrevCluster(clusterUnits[0]);
        }else{
            throw new RuntimeException(""+changeType);
        }



        NewWVTreeLikelihood prevTreeLikelihood = treeLiks.get(prevCluster);
        for(int i = 0; i < clusterUnits.length;i++){
            //Get sites associated to each unit.
            int[] sites = ((GeneralUnitAlignment)alignment).getSitesByUnit(clusterUnits[i]);
            for(int j = 0; j < sites.length; j++){
                int patternIndex = alignment.getPatternIndex(sites[j]);
                patternWeights[patternIndex]++;
                prevTreeLikelihood.removeWeight(patternIndex,1);

            }
        }

        NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(
                patternWeights,
                alignment,
                (Tree) treeInput.get(),
                useAmbiguitiesInput.get(),
                siteModel,
                branchRateModelInput.get());
        try{

            treeLik.calculateLogP();
            treeLik.store();
            treeLiks.add(dpSiteModel.getLastAddedIndex(),treeLik);
        }catch(Exception e){
            throw new RuntimeException(e);
        }



    }

    public void mergeTreeLikelihoods(){

        NewWVTreeLikelihood removedTreeLikelihood =  treeLiks.remove(dpSiteModel.getRemovedIndex());
        NewWVTreeLikelihood mergedTreeLikelihood;
        if(changeType == ChangeType.MERGE || changeType == ChangeType.MERGE_AND_VALUE_CHANGE){
            mergedTreeLikelihood =  treeLiks.get(dpSiteModel.getDirtySiteModelIndex());
        }else if(changeType == ChangeType.REMOVED){
            int dirtyUnit = dpSiteModel.getLastDirtySite();
            mergedTreeLikelihood = treeLiks.get(dpSiteModel.getCurrCluster(dirtyUnit));
        }else{
            throw new RuntimeException();
        }
        int [] patternWeights = removedTreeLikelihood.getPatternWeights();
        for(int i = 0;i < patternWeights.length;i++){
            mergedTreeLikelihood.addWeight(i,patternWeights[i]);
        }
    }

    protected void updateWeights(){
        int dirtyUnit = dpSiteModel.getLastDirtySite();

        if(changeType==ChangeType.POINTER_CHANGED){

            int prevCluster = dpSiteModel.getPrevCluster(dirtyUnit);
            int currCluster = dpSiteModel.getCurrCluster(dirtyUnit);
            int[] sites = ((GeneralUnitAlignment)alignment).getSitesByUnit(dirtyUnit);
            for(int i = 0; i < sites.length;i++){
                treeLiks.get(prevCluster).removeWeight(alignment.getPatternIndex(sites[i]),1);
                treeLiks.get(currCluster).addWeight(alignment.getPatternIndex(sites[i]),1);
            }

        }else if(changeType == ChangeType.POINTERS_SWAPPED){
            int[] swappedUnits = dpSiteModel.getSwappedSites();

            for(int i = 0; i < swappedUnits.length; i++){
                int [] dirtySites = ((GeneralUnitAlignment)alignment).getSitesByUnit(swappedUnits[i]);
                int prevCluster = dpSiteModel.getPrevCluster(swappedUnits[i]);
                int currCluster = dpSiteModel.getCurrCluster(swappedUnits[i]);

                for(int j = 0; j < dirtySites.length; j++){
                    treeLiks.get(prevCluster).removeWeight(alignment.getPatternIndex(dirtySites[j]),1);
                    treeLiks.get(currCluster).addWeight(alignment.getPatternIndex(dirtySites[j]),1);
                }
            }

        }
    }

    private void handlePointersChange(){
        int[] dirtyUnits = dpSiteModel.getLastDirtySites();

        for(int dirtyUnit: dirtyUnits){
            int[] dirtySites = ((GeneralUnitAlignment)alignment).getSitesByUnit(dirtyUnit);
            int prevCluster = dpVal.getPrevCategory(dirtyUnit);
            int currCluster = dpVal.getCurrCategory(dirtyUnit);
            for(int dirtySite: dirtySites){
                treeLiks.get(prevCluster).removeWeight(alignment.getPatternIndex(dirtySite),1);
                treeLiks.get(currCluster).addWeight(alignment.getPatternIndex(dirtySite),1);
            }
        }


    }



    @Override
    protected boolean requiresRecalculation() {
        boolean recalculate = false;
        if(dpSiteModel.isDirtyCalculation()){

            changeType = dpSiteModel.getChangeType();
            //System.out.println("treeLik requires recal!!"+changeType);
            if(changeType == ChangeType.ADDED){
                //System.out.println("added!!");
                splitTreeLikelihood();
                //this.changeType = ChangeType.ADDED;
                updateWeights();

            }else if(changeType == ChangeType.REMOVED){
                //System.out.println("removed!!");
                mergeTreeLikelihoods();
                //this.changeType = ChangeType.REMOVED;

            }else if(changeType == ChangeType.SPLIT){

                splitTreeLikelihood();
                //this.changeType = ChangeType.SPLIT;

            }else if(changeType == ChangeType.SPLIT_AND_VALUE_CHANGE){

                splitTreeLikelihood();
                MCMCNodeFactory.checkDirtiness(treeLiks.get(dpSiteModel.getDirtySiteModelIndex()));
                //this.changeType = ChangeType.SPLIT_AND_VALUE_CHANGE;

            }else if(changeType == ChangeType.MERGE){

                mergeTreeLikelihoods();
                //this.changeType = ChangeType.MERGE;

            }else if(changeType == ChangeType.MERGE_AND_VALUE_CHANGE){

                mergeTreeLikelihoods();
                MCMCNodeFactory.checkDirtiness(treeLiks.get(dpSiteModel.getDirtySiteModelIndex()));
                //this.changeType = ChangeType.MERGE_AND_VALUE_CHANGE;

            }else if (changeType == ChangeType.POINTER_CHANGED||  changeType == ChangeType.POINTERS_SWAPPED){
                //this.changeType = changeType;
                updateWeights();
            }else if(changeType == ChangeType.MULTIPLE_POINTERS_CHANGED){
                handlePointersChange();
                //this.changeType = changeType;
            }else if(changeType == ChangeType.VALUE_CHANGED){
                this.changeType = ChangeType.VALUE_CHANGED;

            }else{
                this.changeType = ChangeType.ALL;
            }



            recalculate = true;
        }else if(treeInput.get().somethingIsDirty()){
            recalculate = true;

        }else if(branchRateModelInput.get().isDirtyCalculation()){
            recalculate = true;
        }
        if(recalculate){

            for(NewWVTreeLikelihood treeLik:treeLiks){

                MCMCNodeFactory.checkDirtiness(treeLik);
            }
        }

        return recalculate;
    }

    public double getSiteLogLikelihood(int iCluster, int unitIndex){

        GeneralUnitAlignment alignment = (GeneralUnitAlignment)this.alignment;
        int[] sites = alignment.getSitesByUnit(unitIndex);
        double logP = 0.0;
        for(int site:sites){
            logP += treeLiks.get(iCluster).getPatternLogLikelihood(alignment.getPatternIndex(site));
        }
        return logP;
    }

}
