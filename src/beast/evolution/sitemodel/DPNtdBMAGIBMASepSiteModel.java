package beast.evolution.sitemodel;

import beast.core.Input;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.ChangeType;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
public class DPNtdBMAGIBMASepSiteModel extends DPNtdRateSepSiteModel {
    public Input<ParameterList> alphaListInput = new Input<ParameterList>(
            "alphaList",
            "A list of unique alpha values of Gamma distribution used to model sites rates.",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> siteModelChoiceListInput = new Input<ParameterList>(
            "siteModelChoiceList",
            "A list of unique indicator values that determines whether the gamma site model includes the alpha or proportion invariant parameter.",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> invPrListInput = new Input<ParameterList>(
        "invPrList",
        "a list of unique invariant proportion values used to model sites rates.",
            Input.Validate.REQUIRED
    );

    public Input<Integer> gammaCategoryCountInput =
            new Input<Integer>("gammaCategoryCount", "gamma category count (default=zero for no gamma)", 1);

    public Input<Boolean> invPrLogitInput = new Input<Boolean>(
            "invPrLogit",
            "Is transforming the invPr to logit space.",
            false
    );

    protected ParameterList alphaList;
    protected ParameterList invPrList;
    protected int gammaCategoryCount;
    protected ParameterList siteModelChoiceList;
    protected boolean invPrLogit;


    public void initAndValidate() throws Exception{
        gammaCategoryCount =  gammaCategoryCountInput.get();
        siteModelChoiceList = siteModelChoiceListInput.get();
        alphaList = alphaListInput.get();
        invPrList = invPrListInput.get();
        invPrLogit = invPrLogitInput.get();
        super.initAndValidate();

    }

    /*
     * Setup the site model and weight matrices.
     */
    public void setup(){
        siteModelsMatrix = new QuietSiteModel[substClusterLimit][ratesClusterLimit];
        storedSiteModelsMatrix = new QuietSiteModel[substClusterLimit][ratesClusterLimit];
        siteModelWeights = new int[substClusterLimit][ratesClusterLimit];
        storedSiteModelWeights = new int[substClusterLimit][ratesClusterLimit];
        //substWeights = new int[substClusterLimit];
        //ratesWeights = new int[ratesClusterLimit];
        clusterMap = new int[2][ratesPointers.getDimension()];
        storedClusterMap = new int[2][ratesPointers.getDimension()];

        int[] substModelPointerIndicies = dpNtdBMA.getPointerIndices();
        int substModelIndex = -1;
        int rateIndex = -1;

        for(int i = 0; i < eltCount; i++){
            substModelIndex = substModelPointerIndicies[i];
            //System.out.println("eltCount: "+substModelIndex);
            rateIndex = ratesPointers.indexInList(i,ratesList);
            if(siteModelsMatrix[substModelIndex][rateIndex] == null){

                try{

                    QuietGammaSiteBMA siteModel = new QuietGammaSiteBMA(
                            dpNtdBMA.getModel(substModelIndex),
                            ratesList.getParameter(rateIndex),
                            alphaList.getParameter(rateIndex),
                            invPrList.getParameter(rateIndex),
                            true,
                            gammaCategoryCount,
                            siteModelChoiceList.getParameter(rateIndex),
                            invPrLogit
                    );
                    siteModelsMatrix[substModelIndex][rateIndex] = siteModel;
                    storedSiteModelsMatrix[substModelIndex][rateIndex] = siteModel;
                    siteModels.add(siteModel);

                }catch(Exception e){
                    throw new RuntimeException(e);
                }

            }
            siteModelWeights[substModelIndex][rateIndex]++;
            storedSiteModelWeights[substModelIndex][rateIndex]++;
            //substWeights[substModelIndex]++;
            //ratesWeights[rateIndex]++;
            clusterMap[NTDBMA][i] = substModelIndex;
            clusterMap[RATES][i] = rateIndex;
            storedClusterMap[NTDBMA][i] = substModelIndex;
            storedClusterMap[RATES][i] = rateIndex;
        }

    }

    /*
     * Handles changes in pointers that does not involve changes in the number of clusters
     */
    public void handlePointerChange(int lastDirtySite) throws Exception{

        //Previous cluster ids of the last dirty site
        int prevNtdBMAIdNum = clusterMap[NTDBMA][lastDirtySite];
        int prevRateIdNum = clusterMap[RATES][lastDirtySite];

        //Current cluster ids of the last dirty site
        SwitchingNtdBMA ntdBMA = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(lastDirtySite)); //todo
        int listIndex = ratesPointers.indexInList(lastDirtySite,ratesList);
        QuietRealParameter muParameter = ratesPointers.getParameter(lastDirtySite);
        QuietRealParameter alphaParameter = alphaList.getParameter(listIndex);
        QuietRealParameter invPrParameter = invPrList.getParameter(listIndex);
        QuietRealParameter siteModelChoice = siteModelChoiceList.getParameter(listIndex);

        //int[] siteClusterMap = new int[2];
        int ntdBMACluster = ntdBMA.getIDNumber();
        int rateCluster = muParameter.getIDNumber();



        //update mapping and weights
        updateMap(lastDirtySite, ntdBMACluster, rateCluster);



        //If the propose combination is new then create a new site model
        if(siteModelsMatrix[ntdBMACluster][rateCluster] == null){
            QuietGammaSiteBMA siteModel = new QuietGammaSiteBMA(
                    ntdBMA,
                    muParameter,
                    alphaParameter,
                    invPrParameter,
                    true,
                    gammaCategoryCount,
                    siteModelChoice,
                    invPrLogit
            );
            siteModelsMatrix[ntdBMACluster][rateCluster] = siteModel;
            siteModels.add(siteModel);
        }

        //If the previous combination no longer has any weight then remove
        //System.out.println(prevNtdBMAIdNum+" "+prevRateIdNum);

    }


    public void updateModelMatrix(int siteIndex) throws Exception{
        //Remove site model if it has zero pattern weight
        if(siteModelWeights[storedClusterMap[NTDBMA][siteIndex]][storedClusterMap[RATES][siteIndex]] == 0
                && siteModelsMatrix[storedClusterMap[NTDBMA][siteIndex]][storedClusterMap[RATES][siteIndex]] != null){
            siteModels.remove(siteModels.indexOf(siteModelsMatrix[storedClusterMap[NTDBMA][siteIndex]][storedClusterMap[RATES][siteIndex]]));
            siteModelsMatrix[storedClusterMap[NTDBMA][siteIndex]][storedClusterMap[RATES][siteIndex]] = null;
        }

        //Add site model if this is a new substModel and rate combination
        if(siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] == null){
            int listIndex = ratesPointers.indexInList(siteIndex,ratesList);
            //System.out.println("add: "+clusterMap[NTDBMA][siteIndex]+" "+clusterMap[RATES][siteIndex]);
            QuietGammaSiteBMA siteModel = new QuietGammaSiteBMA(
                    dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(siteIndex)),
                    ratesPointers.getParameter(siteIndex),
                    alphaList.getParameter(listIndex),
                    invPrList.getParameter(listIndex),
                    true,
                    gammaCategoryCount,
                    siteModelChoiceList.getParameter(listIndex),
                    invPrLogit
            );
            siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] = siteModel;
            siteModels.add(siteModel);
        }

    }


    /*
     * Update cluster mapping
     */
    public void updateMap(int siteIndex, int ntdBMACluster, int rateCluster) throws Exception{
        /*for(int i = 0; i < clusterMap.length;i++){
            System.out.print("NTDBMA: ");
            for(int j = 0; j < clusterMap[i].length;j++){
                System.out.print(clusterMap[i][j]+" ");
            }
            System.out.println();
        }

        for(int i = 0; i < storedClusterMap.length;i++){
            System.out.print("stored NTDBMA: ");
            for(int j = 0; j < storedClusterMap[i].length;j++){
                System.out.print(storedClusterMap[i][j]+" ");
            }
            System.out.println();
        }*/

        mappingChanged = true;

        //Update the weight matrix
        moveWeight(
                clusterMap[NTDBMA][siteIndex],  //current substModel
                clusterMap[RATES][siteIndex],   //current rate
                ntdBMACluster,                  //potentially new substModel
                rateCluster,                    //potentially new rate
                1
        );

        int prevNtdBMAId = clusterMap[NTDBMA][siteIndex];
        int prevRateId = clusterMap[RATES][siteIndex];

        //Update mapping
        clusterMap[NTDBMA][siteIndex] = ntdBMACluster;
        clusterMap[RATES][siteIndex] = rateCluster;

        //Stored map
        //storedClusterMap[NTDBMA][siteIndex] = prevNtdBMAId;
        //storedClusterMap[RATES][siteIndex] = prevRateId;

        //System.out.println("flag1: "+siteModelWeights[prevNtdBMAId][prevRateId]);


        //Remove site model if it has zero pattern weight
        if(siteModelWeights[prevNtdBMAId][prevRateId] == 0){
            //System.out.println("remove");
            //System.out.println("siteModels: "+siteModels.size());
            //System.out.println(siteModelsMatrix[prevNtdBMAId][prevRateId]);
            siteModels.remove(siteModels.indexOf(siteModelsMatrix[prevNtdBMAId][prevRateId]));
            siteModelsMatrix[prevNtdBMAId][prevRateId] = null;
        }

        //Add site model if this is a new substModel and rate combination
        if(siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] == null){
            //System.out.println("add: "+clusterMap[NTDBMA][siteIndex]+" "+clusterMap[RATES][siteIndex]);


            int listIndex = ratesPointers.indexInList(siteIndex,ratesList);
            //System.out.println("add: "+clusterMap[NTDBMA][siteIndex]+" "+clusterMap[RATES][siteIndex]);

            QuietGammaSiteBMA siteModel = new QuietGammaSiteBMA(
                    dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(siteIndex)),
                    ratesPointers.getParameter(siteIndex),
                    alphaList.getParameter(listIndex),
                    invPrList.getParameter(listIndex),
                    true,
                    gammaCategoryCount,
                    siteModelChoiceList.getParameter(listIndex),
                    invPrLogit
            );
            siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] = siteModel;
            siteModels.add(siteModel);
        }

        /*for(int i = 0; i < clusterMap.length;i++){
            System.out.print("RATE: ");
            for(int j = 0; j < clusterMap[i].length;j++){
                System.out.print(clusterMap[i][j]+" ");
            }
            System.out.println();
        }
        for(int i = 0; i < storedClusterMap.length;i++){
            System.out.print("stored RATE: ");
            for(int j = 0; j < storedClusterMap[i].length;j++){
                System.out.print(storedClusterMap[i][j]+" ");
            }
            System.out.println();
        }*/




    }


    public void checkSiteModelsDirtiness(int changedInput){
        int fixedIndex = -1;

        if(changedInput == NTDBMA){

            fixedIndex = dpNtdBMA.getModel(dpNtdBMA.getDirtyModelIndex()).getIDNumber(); //todo use getDirtyModelIDNumber
            for(int i = 0; i < siteModelsMatrix[fixedIndex].length;i++){
                if(siteModelsMatrix[fixedIndex][i] != null){
                    MCMCNodeFactory.checkDirtiness(siteModelsMatrix[fixedIndex][i]);
                    //siteModelsMatrix[fixedIndex][i].checkDirtiness();
                    //System.out.println(i+" "+fixedIndex+" "+siteModelsMatrix[fixedIndex][i].isDirtyCalculation()+" "+siteModelsMatrix[fixedIndex][i].getRateParameter().getIDNumber());
                }
            }
            //System.out.println(getID()+" fixedIndex: "+fixedIndex);
        }else if(changedInput == RATES){
            if(ratesList.somethingIsDirty()){
                fixedIndex = ratesList.getParameter(ratesList.getDirtyIndex()).getIDNumber(); //todo use getDirtyParameterIDNumber
            }else if(alphaList.somethingIsDirty()){
                fixedIndex = alphaList.getDirtyParameterIDNumber();
            }else if(invPrList.somethingIsDirty()){
                fixedIndex = invPrList.getDirtyParameterIDNumber();
            }else{
                fixedIndex = siteModelChoiceList.getDirtyParameterIDNumber();
            }
            for(int i = 0; i < siteModelsMatrix.length;i++){
                if(siteModelsMatrix[i][fixedIndex] != null){
                    MCMCNodeFactory.checkDirtiness(siteModelsMatrix[i][fixedIndex]);
                }
            }

        }else{
            throw new RuntimeException("Can only remove clusters for either ntdBMA model or rates.");
        }

    }


    public boolean requiresRecalculation(){

        boolean recalculate = false;
        mappingChanged = false;
        try{

            //ChangeType substModelChangeType = dpNtdBMA.getChangeType();
            if(ratesList.somethingIsDirty() ||
                    alphaList.somethingIsDirty() ||
                    invPrList.somethingIsDirty() ||
                    siteModelChoiceList.somethingIsDirty() ||
                    ratesPointers.somethingIsDirty()||
                    dpNtdBMA.isDirtyCalculation()){
                //ChangeType inputChangeType = null;


                //Check whether subst model or rates have changed and
                //retrieve change type.
                if(dpNtdBMA.isDirtyCalculation()){

                    changeType = dpNtdBMA.getChangeType();
                    changedInput = NTDBMA;
                    lastDirtySite = dpNtdBMA.getLastDirtySite();

                }else{
                    if(ratesList.somethingIsDirty()){
                        changeType = ratesList.getChangeType();
                    }else if(alphaList.somethingIsDirty()){
                        changeType = alphaList.getChangeType();
                    }else if(invPrList.somethingIsDirty()){
                        changeType = invPrList.getChangeType();
                    }else if(siteModelChoiceList.somethingIsDirty()){
                        changeType = siteModelChoiceList.getChangeType();
                    }else{
                        changeType = ratesPointers.getChangeType();
                    }
                    lastDirtySite = ratesPointers.getLastDirty();
                    changedInput = RATES;
                }


                //System.out.println("changeType: "+this.changeType);
                if(changeType == ChangeType.ADDED){
                    //System.out.println("ADD");
                    addSiteModel(changedInput, lastDirtySite);

                }else if(changeType == ChangeType.SPLIT){

                    handleSplit(changedInput);

                }else if(changeType == ChangeType.REMOVED){
                    //System.out.println("REMOVED");
                    removeSiteModel(changedInput, lastDirtySite);

                }else if(changeType == ChangeType.MERGE){

                    handleMerge(changedInput);

                }else if(changeType == ChangeType.POINTER_CHANGED){
                    //System.out.println("POINTER_CHANGED");

                    handlePointerChange(lastDirtySite);
                }else if(changeType == ChangeType.MULTIPLE_POINTERS_CHANGED){

                    handleMultiPointerChanges(changedInput);


                }else if(changeType == ChangeType.POINTERS_SWAPPED){
                    if(changedInput == NTDBMA){
                        clusterSites = dpNtdBMA.getSwappedSites();

                    }else{
                        clusterSites = ratesPointers.getSwappedSites();
                    }

                    if(clusterSites[0] != clusterSites[1]){
                        handlePointerSwap(clusterSites[0], clusterSites[1]);
                        //handlePointerChange(clusterSites[0]);
                        //handlePointerChange(clusterSites[1]);

                    }

                }else if(changeType == ChangeType.VALUE_CHANGED){
                    //System.out.println("VALUE_CHANGED");
                    //When there's a value change, the cluster assignment doesn't change.
                    checkSiteModelsDirtiness(changedInput);

                }else {
                    this.changeType = ChangeType.ALL;
                    for(int i = 0; i < siteModelsMatrix.length; i++){
                        for(int j = 0; j < siteModelsMatrix[i].length; j++){
                            if(siteModelsMatrix[i][j] != null){
                                MCMCNodeFactory.checkDirtiness(siteModelsMatrix[i][j]);

                            }

                        }
                    }
                    //setupPointerIndices();
                }
                recalculate = true;
                //System.err.println("dirty1");

            }



        }catch(Exception e){
            throw new RuntimeException(e);
        }
        //System.out.println("siteModel changeType: "+changeType);
        return recalculate;
    }

    public double getRateValue(int unitIndex){
        return ((QuietGammaSiteBMA)siteModelsMatrix[clusterMap[NTDBMA][unitIndex]][clusterMap[RATES][unitIndex]]).getMuValue();
    }

    public double getAlphaValue(int unitIndex){
        return ((QuietGammaSiteBMA)siteModelsMatrix[clusterMap[NTDBMA][unitIndex]][clusterMap[RATES][unitIndex]]).getAlphaValue();
    }

    public double getInvPrValue(int unitIndex){
        return ((QuietGammaSiteBMA)siteModelsMatrix[clusterMap[NTDBMA][unitIndex]][clusterMap[RATES][unitIndex]]).getInvPrValue();
    }

    public double getModelChoiceValue(int unitIndex){
        return ((QuietGammaSiteBMA)siteModelsMatrix[clusterMap[NTDBMA][unitIndex]][clusterMap[RATES][unitIndex]]).getModelChoiceValue();
    }




    public int getRateID(int unitIndex){


        return siteModelsMatrix[clusterMap[NTDBMA][unitIndex]][clusterMap[RATES][unitIndex]].getRateParameter().getIDNumber();
    }

    public void store(){
        for(SiteModel siteModel:siteModels){
            siteModel.store();

        }
        super.store();
    }

    public void restore(){
        //System.out.println("restore");
        super.restore();
        for(SiteModel siteModel:siteModels){
            siteModel.restore();
        }
    }


}
