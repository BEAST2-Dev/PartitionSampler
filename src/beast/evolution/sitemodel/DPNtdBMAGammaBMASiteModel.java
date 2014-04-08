package beast.evolution.sitemodel;

import beast.core.Input;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.ChangeType;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.evolution.substitutionmodel.DPNtdBMA;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
public class DPNtdBMAGammaBMASiteModel extends DPNtdRateSiteModel{
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

    private ParameterList alphaList;
    private ParameterList invPrList;
    private int gammaCategoryCount;
    private ParameterList siteModelChoiceList;
    private boolean invPrLogit;
    public void initAndValidate() throws Exception{
        //super.initAndValidate();
        dpNtdBMA = dpNtdBMAInput.get();

        int ntdBMACount = dpNtdBMA.getDimension();

        ratePointers = ratePointersInput.get();
        rateList = rateListInput.get();
        alphaList = alphaListInput.get();
        invPrList = invPrListInput.get();

        if(rateList.getDimension() != dpNtdBMA.getDimension()){
            throw new RuntimeException("The number of clustres for rates and substution models must be the same.");
        }

        //Setting up site models
        siteModels = new ArrayList<QuietSiteModel>();
        gammaCategoryCount =  gammaCategoryCountInput.get();
        System.out.println(getID()+" "+gammaCategoryCount);

        siteModelChoiceList = siteModelChoiceListInput.get();
        invPrLogit = invPrLogitInput.get();
        for(int i = 0;i < ntdBMACount; i++){
            QuietRealParameter muParameter = rateList.getParameter(ratePointers.indexInList(i,rateList));
            QuietRealParameter shapeParameter = alphaList.getParameter(i);
            QuietRealParameter invPrParameter = invPrList.getParameter(i);
            QuietRealParameter modelChoice = siteModelChoiceList.getParameter(i);
            QuietSiteModel siteModel = new QuietGammaSiteBMA(
                    dpNtdBMA.getModel(i),
                    muParameter,
                    shapeParameter,
                    invPrParameter,
                    true,
                    gammaCategoryCount,
                    modelChoice,
                    invPrLogit
            );
            siteModels.add(siteModel);
        }
    }


    protected void addSiteModel(){
        try{

            QuietRealParameter muParameter = rateList.getParameter(rateList.getLastAddedIndex());
            QuietRealParameter shapeParameter = alphaList.getParameter(alphaList.getLastAddedIndex());
            QuietRealParameter invPrParameter = invPrList.getParameter(invPrList.getLastAddedIndex());
            QuietRealParameter modelChoice = siteModelChoiceList.getParameter(siteModelChoiceList.getLastAddedIndex());
            QuietGammaSiteBMA siteModel = new QuietGammaSiteBMA(
                    dpNtdBMA.getModel(dpNtdBMA.getLastAddedIndex()),
                    muParameter,
                    shapeParameter,
                    invPrParameter,
                    true,
                    gammaCategoryCount,
                    modelChoice,
                    invPrLogit
                    );
            siteModels.add(dpNtdBMA.getLastAddedIndex(),siteModel);
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }


    public boolean requiresRecalculation(){

        boolean recalculate = false;
        //System.err.println("dirty0");
        ChangeType substModelChangeType = dpNtdBMA.getChangeType();
        //System.out.println(rateList.somethingIsDirty() +" "+ dpNtdBMA.isDirtyCalculation());
        //System.out.println(changeType);
        if(rateList.somethingIsDirty() && dpNtdBMA.isDirtyCalculation()){
            changeType = rateList.getChangeType();

            if(changeType != substModelChangeType){
                System.out.println(substModelChangeType+" "+changeType);
                throw new RuntimeException("Can only handle same type of changes to subst and rate at one time.");
            }


            if(changeType == ChangeType.ADDED||changeType == ChangeType.SPLIT){

                addSiteModel();

            }else if(changeType == ChangeType.REMOVED || changeType == ChangeType.MERGE){

                removeSiteModel(rateList.getRemovedIndex());

            }else if(changeType == ChangeType.VALUE_CHANGED){

                for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }

            }else{

                this.changeType = ChangeType.ALL;
                for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }

            }
            recalculate = true;

        }else if (ratePointers.somethingIsDirty()){

            changeType = ratePointers.getChangeType();
            if(changeType != substModelChangeType){
                System.out.println(changeType+" "+substModelChangeType);
                throw new RuntimeException("Can only handle same type of changes to subst and rate at one time.");
            }
            recalculate = true;

        }else{
            if(rateList.somethingIsDirty()){
                this.changeType = rateList.getChangeType();
                recalculate = true;

            }else if(alphaList.somethingIsDirty()){
                this.changeType = alphaList.getChangeType();
                recalculate = true;

            }else if(invPrList.somethingIsDirty()){
                this.changeType = invPrList.getChangeType();
                recalculate = true;

            }else if(siteModelChoiceList.somethingIsDirty()){
                this.changeType = siteModelChoiceList.getChangeType();
                recalculate = true;

            }else if(dpNtdBMA.isDirtyCalculation() && dpNtdBMA.getChangeType() == ChangeType.VALUE_CHANGED){
                this.changeType = ChangeType.VALUE_CHANGED;

                recalculate = true;
            }
            for(SiteModel siteModel:siteModels){

                MCMCNodeFactory.checkDirtiness(siteModel);

            }
        }

        return recalculate;
    }


    public void store(){
        //System.out.println("store");

        for(QuietSiteModel siteModel:siteModels){
            siteModel.store();
        }
        super.store();
    }

    public void restore(){
        //System.out.println("restore");
        super.restore();
        for(QuietSiteModel siteModel:siteModels){
            siteModel.restore();
        }
    }

}
