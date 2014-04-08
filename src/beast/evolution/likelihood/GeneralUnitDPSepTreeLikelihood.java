package beast.evolution.likelihood;

import beast.app.BeastMCMC;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.ChangeType;
import beast.evolution.alignment.GeneralUnitAlignment;
import beast.evolution.sitemodel.DPNtdRateSepSiteModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 15/07/13
 * Time: 3:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneralUnitDPSepTreeLikelihood extends DPSepTreeLikelihood {

    public void initAndValidate() throws Exception{
        useThreads = useThreadsInput.get() && (BeastMCMC.m_nThreads > 1);
        useThreadsEvenly = useThreadsEvenlyInput.get() && (BeastMCMC.m_nThreads > 1);

        alignment = m_data.get();
        int patternCount = alignment.getPatternCount();
        if(!(m_pSiteModel.get() instanceof DPNtdRateSepSiteModel)){
            throw new RuntimeException("DPNtdRateSepSiteModel required for site model.");
        }
        dpSiteModel = (DPNtdRateSepSiteModel)m_pSiteModel.get();
        int siteModelCount = dpSiteModel.getSiteModelCount();

        /*
         * Set up a 3 dimensional array to store the weights of each ntdBMA/rate combination (dimensions: ntdBMA, rates, weights).
         */
        int[][][]clusterWeights = new int[dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()][patternCount];
        int[] clusterIds;
        int siteCount = alignment.getSiteCount();
        for(int i = 0; i < siteCount; i++){
            clusterIds = dpSiteModel.getCurrClusters(((GeneralUnitAlignment)alignment).getUnitBySite(i));
            clusterWeights[clusterIds[DPNtdRateSepSiteModel.NTDBMA]][clusterIds[DPNtdRateSepSiteModel.RATES]][alignment.getPatternIndex(i)]++;
        }

        /*
         *  Traverse through the list site models and create tree likelihoods.
         *  The tree likelihoods are then added to a list and the matrix.
         *  The order of likelihoods in the list corresponds to that of site models in the DPSiteModel list.
         */
        //treeLiksMatrix = new WVTreeLikelihood[dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()];
        treeLiksMatrix = new NewWVTreeLikelihood[dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()];
        storedTreeLiksMatrix = new NewWVTreeLikelihood[dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()];
        for(int i = 0; i < siteModelCount;i++){

            //Get the ids and hence the positions of siteModel/treeLikelihoods in a list
            int ntdBMAId = ((SwitchingNtdBMA)dpSiteModel.getSiteModel(i).getSubstitutionModel()).getIDNumber();
            int ratesId = dpSiteModel.getSiteModel(i).getRateParameter().getIDNumber();

            //Create the tree likelihood
            //WVTreeLikelihood treeLik = new WVTreeLikelihood(clusterWeights[ntdBMAId][ratesId]);

            NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(
                    clusterWeights[ntdBMAId][ratesId],
                    alignment,
                    m_tree.get(),
                    useAmbiguitiesInput.get(),
                    dpSiteModel.getSiteModel(i),
                    m_pBranchRateModel.get());

            //Add to list and matrix for the convenience of processesing
            treeLiks.add(treeLik);
            treeLiksMatrix[ntdBMAId][ratesId] = treeLik;
            storedTreeLiksMatrix[ntdBMAId][ratesId] = treeLik;
        }
        if(!(alignment instanceof GeneralUnitAlignment)){
            throw new RuntimeException("Must use GeneralUnitAlignment object.");
        }
    }


    protected void update(){
        int dirtyUnit = dpSiteModel.getLastDirtySite();
        int[] currCategoryIds = dpSiteModel.getCurrClusters(dirtyUnit);
        int[] prevCategoryIds = dpSiteModel.getPrevClusters(dirtyUnit);

        //Create a new tree likelihood (with zero weights) when there is a new site model.
        if(treeLiksMatrix[currCategoryIds[DPNtdRateSepSiteModel.NTDBMA]][currCategoryIds[DPNtdRateSepSiteModel.RATES]] == null){
            addTreeLikelihood(currCategoryIds[DPNtdRateSepSiteModel.NTDBMA],currCategoryIds[DPNtdRateSepSiteModel.RATES]);
        }

        //Move weight

        int[] dirtySites = ((GeneralUnitAlignment)alignment).getSitesByUnit(dirtyUnit);
        for(int dirtySite:dirtySites){
            moveWeight(
                prevCategoryIds[DPNtdRateSepSiteModel.NTDBMA],
                prevCategoryIds[DPNtdRateSepSiteModel.RATES],
                currCategoryIds[DPNtdRateSepSiteModel.NTDBMA],
                currCategoryIds[DPNtdRateSepSiteModel.RATES],
                dirtySite,
                1
            );
        }

        //Remove likelihoods that have zero weights
        if(dpSiteModel.getCombinationWeight(prevCategoryIds[DPNtdRateSepSiteModel.NTDBMA], prevCategoryIds[DPNtdRateSepSiteModel.RATES]) == 0){
            treeLiks.remove(treeLiksMatrix[prevCategoryIds[DPNtdRateSepSiteModel.NTDBMA]][prevCategoryIds[DPNtdRateSepSiteModel.RATES]]);
            treeLiksMatrix[prevCategoryIds[DPNtdRateSepSiteModel.NTDBMA]][prevCategoryIds[DPNtdRateSepSiteModel.RATES]] = null;
        }
    }


    protected void updates(){
        int[] dirtyUnits =  dpSiteModel.getLastDirtySites();
        for(int dirtyUnit:dirtyUnits){
            update(dirtyUnit);
        }

        for(int dirtyUnit:dirtyUnits){
            int[] prevClusterIds = dpSiteModel.getPrevClusters(dirtyUnit);

            //Remove likelihoods that have zero weights
            if(dpSiteModel.getCombinationWeight(prevClusterIds[DPNtdRateSepSiteModel.NTDBMA], prevClusterIds[DPNtdRateSepSiteModel.RATES]) == 0){
                treeLiks.remove(treeLiksMatrix[prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]]);
                treeLiksMatrix[prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]] = null;
            }
        }
    }

    protected void update(int dirtyUnit){

        int[] currCategoryIds = dpSiteModel.getCurrClusters(dirtyUnit);
        int[] prevCategoryIds = dpSiteModel.getPrevClusters(dirtyUnit);

        //Create a new tree likelihood (with zero weights) when there is a new site model.
        if(treeLiksMatrix[currCategoryIds[DPNtdRateSepSiteModel.NTDBMA]][currCategoryIds[DPNtdRateSepSiteModel.RATES]] == null){
            addTreeLikelihood(currCategoryIds[DPNtdRateSepSiteModel.NTDBMA],currCategoryIds[DPNtdRateSepSiteModel.RATES]);
        }


        //Get sites associated with the dirty unit
        int[] dirtySites = ((GeneralUnitAlignment)alignment).getSitesByUnit(dirtyUnit);
        for(int dirtySite:dirtySites){
            //Move weight
            moveWeight(
                prevCategoryIds[DPNtdRateSepSiteModel.NTDBMA],
                prevCategoryIds[DPNtdRateSepSiteModel.RATES],
                currCategoryIds[DPNtdRateSepSiteModel.NTDBMA],
                currCategoryIds[DPNtdRateSepSiteModel.RATES],
                dirtySite,
                1
            );
        }


    }

    public void addTreeLikelihood(int substModelID, int rateID){
        //SiteModel siteModel = dpSiteModel.getLastAdded();
        SiteModel siteModel = dpSiteModel.getSiteModel(substModelID, rateID);


        int[] patternWeights = new int[alignment.getPatternCount()];

        //WVTreeLikelihood treeLik = new WVTreeLikelihood(patternWeights);
        //NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(patternWeights);
        NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(
                    patternWeights,
                    alignment,
                    m_tree.get(),
                    useAmbiguitiesInput.get(),
                    siteModel,
                    m_pBranchRateModel.get());
        try{



            treeLik.calculateLogP();
            treeLik.store();
            treeLiks.add(treeLik);
            treeLiksMatrix[substModelID][rateID] = treeLik;
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }


    @Override
    protected boolean requiresRecalculation() {



        boolean recalculate = false;
        if(dpSiteModel.isDirtyCalculation()){

            changeType = dpSiteModel.getChangeType();
            //System.out.println("treeLik requires recal!!"+changeType);
            if(changeType == ChangeType.ADDED || changeType == ChangeType.REMOVED || changeType == ChangeType.POINTER_CHANGED){
                //System.out.println("changeType: "+changeType);
                update();
            }else if(changeType == ChangeType.SPLIT || changeType == ChangeType.MERGE || changeType == ChangeType.MULTIPLE_POINTERS_CHANGED){
                //storeTreeLikelihoods();
                //System.out.println(changeType);
                updates();
            }else if(changeType == ChangeType.POINTERS_SWAPPED){
                int[] dirtySites = dpSiteModel.getLastDirtySites();

                if(dirtySites[0] !=  dirtySites[1]){

                    updates();
                }
            }else if(changeType == ChangeType.VALUE_CHANGED){

                changeType = ChangeType.VALUE_CHANGED;

            }else{
                changeType = ChangeType.ALL;
            }
            recalculate = true;
        }else if(m_tree.get().somethingIsDirty()){
            recalculate = true;

        }else if(m_pBranchRateModel.get().isDirtyCalculation()){
            recalculate = true;
        }

        if(recalculate){
            //SiteModel siteModel = (SiteModel)treeLiks.get(0).m_siteModel;
            for(TreeLikelihood treeLik:treeLiks){
                MCMCNodeFactory.checkDirtiness(treeLik);
                //System.out.println(treeLik.m_siteModel.isDirtyCalculation());

            }


        }

        /*q = dpSiteModel.getSiteModel(0,2);
        if(q != null && treeLiksMatrix[0][2] !=null){

              if(treeLiksMatrix[0][2].getSiteModel() != q){
                  System.out.println("whoa!");
                  treeLiksMatrix[0][2].printThings();
                  throw new RuntimeException("whoa!");
              }


        }    */
        return recalculate;
    }

    public double getSiteLogLikelihood(int inputType, int clusterID, int unitIndex){
        //System.out.println(getClass()+": "+unitIndex);
        GeneralUnitAlignment alignment = (GeneralUnitAlignment)this.alignment;

        double logP = 0.0;
        //System.out.println("pattern: "+alignment.getPatternIndex(siteIndex));
        int[] currClusters = dpSiteModel.getCurrClusters(unitIndex);
        //System.out.println("clusterID: "+clusterID+" "+prevClusters[DPNtdRateSepSiteModel.RATES]+" "+alignment.getPatternIndex(siteIndex));
        if(inputType == DPNtdRateSepSiteModel.NTDBMA){
            if(treeLiksMatrix[clusterID][currClusters[DPNtdRateSepSiteModel.RATES]] != null){

                NewWVTreeLikelihood tmpTL = treeLiksMatrix[clusterID][currClusters[DPNtdRateSepSiteModel.RATES]];
                int[] sites = alignment.getSitesByUnit(unitIndex);
                for(int site:sites){
                    //System.out.println("site: "+site+", "+alignment.getPatternIndex(site)+", "+tmpTL.getPatternLogLikelihood(alignment.getPatternIndex(site)));
                    logP += tmpTL.getPatternLogLikelihood(alignment.getPatternIndex(site));
                }

                return logP;

            }
        }else{

            if(treeLiksMatrix[currClusters[DPNtdRateSepSiteModel.NTDBMA]][clusterID] != null){
                //WVTreeLikelihood tmpTL = treeLiksMatrix[prevClusters[DPNtdRateSepSiteModel.NTDBMA]][clusterID];
                //System.out.println("hi!!");
                NewWVTreeLikelihood tmpTL = treeLiksMatrix[currClusters[DPNtdRateSepSiteModel.NTDBMA]][clusterID];
                int[] sites = alignment.getSitesByUnit(unitIndex);
                for(int site:sites){
                    logP += tmpTL.getPatternLogLikelihood(alignment.getPatternIndex(site));
                }
                return logP;
            }
        }

        return Double.NaN;

    }
}
