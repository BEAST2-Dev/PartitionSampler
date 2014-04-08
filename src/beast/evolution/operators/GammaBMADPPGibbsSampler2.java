package beast.evolution.operators;

import beast.core.Input;
import beast.core.Loggable;
import beast.core.Operator;
import beast.core.parameter.*;
import beast.evolution.likelihood.DPSepTreeLikelihood;
import beast.evolution.likelihood.GeneralUnitSepTempTreeLikelihood;
import beast.evolution.sitemodel.DPNtdRateSepSiteModel;
import beast.math.distributions.CompoundDirichletProcess;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

import java.io.PrintStream;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 10/09/13
 * Time: 6:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class GammaBMADPPGibbsSampler2 extends Operator implements Loggable {


    public Input<DPPointer> ratesPointersInput = new Input<DPPointer>(
            "ratesPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> ratesListInput = new Input<ParameterList>(
            "ratesList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> alphaListInput = new Input<ParameterList>(
            "alphaList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> invPrListInput = new Input<ParameterList>(
            "invPrList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> siteModelListInput = new Input<ParameterList>(
            "siteModelList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );


    public Input<CompoundDirichletProcess> dpInput = new Input<CompoundDirichletProcess>(
            "dirichletProcess",
            "An object of Dirichlet Process",
            Input.Validate.REQUIRED
    );




    public Input<ParametricDistribution> ratesBaseDistrInput = new Input<ParametricDistribution>(
            "ratesBaseDistr",
            "The base distribution of the overall subsitution rate of a site.",
            Input.Validate.REQUIRED
    );


    public Input<ParametricDistribution> alphaBaseDistrInput = new Input<ParametricDistribution>(
            "alphaBaseDistr",
            "The base distribution of the shape parameters of gamma site models",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> invPrBaseDistrInput = new Input<ParametricDistribution>(
            "invPrBaseDistr",
            "The base distribution of the invariant proportion of gamma site models",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> siteModelBaseDistrInput = new Input<ParametricDistribution>(
            "siteModelBaseDistr",
            "The base distribution of the gamma site model",
            Input.Validate.REQUIRED
    );


    public Input<Integer> sampleSizeInput = new Input<Integer>(
            "sampleSize",
            "The number of prelimiary proposals",
            Input.Validate.REQUIRED
    );


    public Input<GeneralUnitSepTempTreeLikelihood> tempLikelihoodInput = new Input<GeneralUnitSepTempTreeLikelihood>(
            "tempLikelihood",
            "The temporary likelihood given the data at site i",
            Input.Validate.REQUIRED
    );

    public Input<DPValuable> dpValuableInput = new Input<DPValuable>(
            "dpVal",
            "reports the counts in each cluster",
            Input.Validate.REQUIRED
    );

    public Input<DPSepTreeLikelihood> dpTreeLikelihoodInput = new Input<DPSepTreeLikelihood>(
            "dpTreeLik",
            "Tree likelihood that handle DPP",
            Input.Validate.REQUIRED
    );

    public Input<Boolean> testCorrectInput = new Input<Boolean>(
            "testCorrect",
            "Whether to check the likelihood calculations are consistent.",
            false
    );

    int tempCount = 1;
    private CompoundDirichletProcess dp;
    private int sampleSize;
    private ParametricDistribution alphaBaseDistr;
    private ParametricDistribution invPrBaseDistr;
    private ParametricDistribution siteModelBaseDistr;
    private ParametricDistribution ratesBaseDistr;
    private DPValuable dpVal;
    private DPSepTreeLikelihood dpTreeLikelihood;
    private int addClusterCount =0;
    private boolean testCorrect;
    public void initAndValidate(){
        testCorrect = testCorrectInput.get();
        dp = dpInput.get();
        ratesBaseDistr = ratesBaseDistrInput.get();
        alphaBaseDistr = alphaBaseDistrInput.get();
        invPrBaseDistr = invPrBaseDistrInput.get();
        siteModelBaseDistr = siteModelBaseDistrInput.get();


        sampleSize = sampleSizeInput.get();
        dpVal = dpValuableInput.get();
        dpTreeLikelihood = dpTreeLikelihoodInput.get();

    }


    public double proposal(){
        //Get the pointer and the list of unique values

        DPPointer ratesPointers = ratesPointersInput.get();
        ParameterList ratesList = ratesListInput.get();
        ParameterList alphaList = alphaListInput.get();
        ParameterList invPrList = invPrListInput.get();
        ParameterList siteModelList = siteModelListInput.get();


        //Randomly pick an index to update, gets it's current value and its position in the parameter list
        int dimPointer = ratesPointers.getDimension();
        int index = Randomizer.nextInt(dimPointer);



        int listIndex = ratesPointers.indexInList(index,ratesList);

        GeneralUnitSepTempTreeLikelihood tempLik = tempLikelihoodInput.get();

        //Get the dimension of the parameter
        int ratesDimValue = ratesList.getParameterDimension();
        RealParameter curr = ratesPointers.getParameter(index);

        //Count the number of items in each cluster but excluding the one about to be updated
        int[] clusterCounts = dpVal.getClusterCounts();
        clusterCounts[listIndex] =  clusterCounts[listIndex]-1;

        QuietRealParameter[] existingRatesVals = new QuietRealParameter[clusterCounts.length];
        QuietRealParameter[] existingAlphaVals = new QuietRealParameter[clusterCounts.length];
        QuietRealParameter[] existingInvPrVals = new QuietRealParameter[clusterCounts.length];
        QuietRealParameter[] existingSiteModel = new QuietRealParameter[clusterCounts.length];

        int[] existingCluster = new int[clusterCounts.length];

        int counter = 0;
        int zeroCount = -1;
        //System.out.println("clusterCounts.length: "+dpVal.getDimension());
        for(int i = 0; i < clusterCounts.length;i++){
            if(clusterCounts[i]>0){
                clusterCounts[counter] = clusterCounts[i];
                existingRatesVals[counter] = ratesList.getParameter(i);
                existingAlphaVals[counter] = alphaList.getParameter(i);
                existingInvPrVals[counter] = invPrList.getParameter(i);
                existingSiteModel[counter] = siteModelList.getParameter(i);
                existingCluster[counter] = i;
                counter++;

            }else{
                zeroCount = i;
            }
        }

        try{

            //Generate a sample of proposals
            QuietRealParameter[] ratesPreProposals = getSamples(ratesBaseDistr, ratesList.getUpper(), ratesList.getLower());
            QuietRealParameter[] alphaPreProposals = getSamples(alphaBaseDistr, alphaList.getUpper(), alphaList.getLower());
            QuietRealParameter[] invPrPreProposals = getSamples(invPrBaseDistr, invPrList.getUpper(), invPrList.getLower());
            QuietRealParameter[] siteModelPreProposals = getSamples(siteModelBaseDistr, siteModelList.getUpper(), siteModelList.getLower());

            //System.err.println("zero count: "+zeroCount);
            //If the a singleton has been picked
            if(zeroCount > -1){
                ratesPreProposals[0] = ratesList.getParameter(zeroCount);
                alphaPreProposals[0] = alphaList.getParameter(zeroCount);
                invPrPreProposals[0] = invPrList.getParameter(zeroCount);
                siteModelPreProposals[0] = siteModelList.getParameter(zeroCount);

            }

            //int dimList = paramList.getDimension();
            int i;
            double concVal =dp.getConcParameter();

            double[] logFullCond = new double[counter+sampleSize];

            for(i = 0; i < counter; i++){

                logFullCond[i] = Math.log(clusterCounts[i]/(dimPointer - 1 + concVal));
                //System.out.println("Complete calculation: ");

                if(testCorrect){
                    double temp1 = tempLik.calculateLogP(
                            existingAlphaVals[i].getValue(),
                            existingInvPrVals[i].getValue(),
                            existingRatesVals[i].getValue(),
                            existingSiteModel[i].getValue(),
                            index
                    );
                    double temp2 = getSiteLogLikelihood(
                            existingAlphaVals[i].getValue(),
                            existingInvPrVals[i].getValue(),
                            existingRatesVals[i].getValue(),
                            existingSiteModel[i].getValue(),
                            ratesList.getParameterIDNumber(existingCluster[i]),
                            index,
                            tempLik
                    );
                    if(temp1 != temp2 && !Double.isNaN(temp2)){
                        System.out.println("existingAlphaVals: "+existingAlphaVals[i]);
                        System.out.println("existingInvPrVals: "+existingInvPrVals[i]);
                        System.out.println("existingRatesVals: "+existingRatesVals[i]);
                        System.out.println("existingSiteModel: "+existingSiteModel[i]);
                        System.out.println(temp1+" "+temp2);
                        throw new RuntimeException();
                    }
                }
                //System.err.println("clusterCounts: "+clusterCounts[i]);
                //logFullCond[i] = logFullCond[i]+ dpTreeLikelihood.getSiteLogLikelihood(existingCluster[i],index);
                logFullCond[i] = logFullCond[i]+getSiteLogLikelihood(
                        existingAlphaVals[i].getValue(),
                        existingInvPrVals[i].getValue(),
                        existingRatesVals[i].getValue(),
                        existingSiteModel[i].getValue(),
                        ratesList.getParameterIDNumber(existingCluster[i]),
                        index,
                        tempLik
                );
            }

            if(zeroCount > -1){
                logFullCond[i] = Math.log(concVal/sampleSize/(dimPointer - 1 + concVal));
                if(testCorrect){
                    double temp1 = tempLik.calculateLogP(
                                    existingAlphaVals[0].getValue(),
                                    existingInvPrVals[0].getValue(),
                                    existingRatesVals[0].getValue(),
                                    existingSiteModel[0].getValue(),
                                    index
                            );
                    //double temp2 = dpTreeLikelihood.getSiteLogLikelihood(zeroCount,index);
                    double temp2 = getSiteLogLikelihood(
                            existingAlphaVals[0].getValue(),
                            existingInvPrVals[0].getValue(),
                            existingRatesVals[0].getValue(),
                            existingSiteModel[0].getValue(),
                            ratesList.getParameterIDNumber(existingCluster[0]),
                            index,
                            tempLik
                    );
                    if(temp1 != temp2){
                        throw new RuntimeException(temp1+" "+temp2);
                    }
                }
                //logFullCond[i] = logFullCond[i]+ dpTreeLikelihood.getSiteLogLikelihood(zeroCount,index);
                logFullCond[i] = logFullCond[i]+getSiteLogLikelihood(
                        alphaPreProposals[0].getValue(),
                        invPrPreProposals[0].getValue(),
                        ratesPreProposals[0].getValue(),
                        siteModelPreProposals[0].getValue(),
                        ratesList.getParameterIDNumber(zeroCount),
                        index,
                        tempLik
                );
                i++;
            }
            for(; i < logFullCond.length; i++){

                logFullCond[i] = Math.log(concVal/sampleSize/(dimPointer - 1 + concVal));
                /*double temp1 = tempLik.calculateLogP(
                                paramPreProposals[i-counter],
                                modelPreProposals[i-counter],
                                freqsPreProposals[i-counter],
                                index
                        );*/

                logFullCond[i] = logFullCond[i]+tempLik.calculateLogP(
                        alphaPreProposals[i-counter].getValue(),
                        invPrPreProposals[i-counter].getValue(),
                        ratesPreProposals[i-counter].getValue(),
                        siteModelPreProposals[i-counter].getValue(),
                        index
                );

                if(Double.isNaN(logFullCond[i])){

                    System.err.println(logFullCond[i]+" "+logFullCond[i]);
                    dpTreeLikelihoodInput.get().getTree().log(-1,System.err);
                    return Double.NEGATIVE_INFINITY;

                }

            }




            //System.err.println("smallestVal2: "+smallestVal);

            double[] fullConditional = new double[logFullCond.length];

            for(i = 0; i < fullConditional.length;i++){
                fullConditional[i] = Math.exp(logFullCond[i]);

                //System.out.println("fullConditional[i]: "+fullConditional[i]+" logFullCond[i]: "+logFullCond[i]);
                if(Double.isNaN(fullConditional[i])){
                    return Double.NEGATIVE_INFINITY;
                }else if(fullConditional[i] == Double.NEGATIVE_INFINITY){

                    fullConditional[i] = 0.0;

                    //return Double.NEGATIVE_INFINITY;
                    //throw new RuntimeException("Huh?");
                    //System.out.println("reject directly, existing parameter vals: "+existingParamVals.length);
                }else if(fullConditional[i] == Double.POSITIVE_INFINITY){
                    //System.out.println("accept directly, existing parameter vals: "+existingParamVals.length);
                    fullConditional[i] = 1.0;

                }
            }

            double total = Randomizer.getTotal(fullConditional);
            if(total < 1e-100){
                /*System.out.println("flag1");
                for(int j = 0; j < fullConditional.length; j++){
                    //fullConditional[j] = Math.exp(logFullCond[j] - smallestVal);
                    System.out.println(fullConditional[j]+" "+logFullCond[j]);
                }*/
                double maxVal = logFullCond[0];
                for(int j = 1; j < logFullCond.length;j++){
                    if(maxVal < logFullCond[j]) {
                        maxVal = logFullCond[j];
                    }
                }
                for(int j = 0; j < fullConditional.length; j++){
                    fullConditional[j] = Math.exp(logFullCond[j] - maxVal);
                    //System.out.println(logFullCond[j] - smallestVal);
                }

            }

            /*for(int j = 0; j < fullConditional.length;j++){
                System.out.println(j+" "+fullConditional[j]);
            }*/


            int proposedIndex = Randomizer.randomChoicePDF(fullConditional);
            //System.out.println("proposedIndex: "+proposedIndex);
            QuietRealParameter ratesProposal;
            QuietRealParameter alphaProposal;
            QuietRealParameter invPrProposal;
            QuietRealParameter siteModelProposal;

            if(proposedIndex < counter){
                alphaProposal = existingAlphaVals[proposedIndex];
                invPrProposal = existingInvPrVals[proposedIndex];
                ratesProposal = existingRatesVals[proposedIndex];
                siteModelProposal = existingSiteModel[proposedIndex];

            }else{
                alphaProposal = alphaPreProposals[proposedIndex-counter];
                invPrProposal = invPrPreProposals[proposedIndex-counter];
                ratesProposal = ratesPreProposals[proposedIndex-counter];
                siteModelProposal = siteModelPreProposals[proposedIndex-counter];
            }

            if(curr != ratesProposal){//Takes a different value from the current
                //System.out.println("before sampler: "+ratesList.getID()+" "+ ratesList.getDimension());
                if(proposedIndex >= counter && zeroCount > -1){//Singleton takes new value

                    alphaList = alphaListInput.get(this);
                    invPrList = invPrListInput.get(this);
                    ratesList = ratesListInput.get(this);
                    siteModelList = siteModelListInput.get(this);

                    int ratesListIndex = ratesPointers.indexInList(index,ratesList);
                    alphaList.setValue(ratesListIndex,0,alphaProposal.getValue());
                    invPrList.setValue(ratesListIndex,0,invPrProposal.getValue());
                    ratesList.setValue(ratesListIndex,0,ratesProposal.getValue());
                    siteModelList.setValue(ratesListIndex,0,siteModelProposal.getValue());
                    zeroCount = -1;

                }else{
                    //Singleton takes existing value or
                    //non-singleton takes new or existing value
                    ratesPointers = ratesPointersInput.get(this);
                    ratesPointers.point(index, ratesProposal);

                    //Non singleton takes new value
                    if(proposedIndex >= counter){
                        //addClusterCount++;
                        alphaList = alphaListInput.get(this);
                        invPrList = invPrListInput.get(this);
                        ratesList = ratesListInput.get(this);
                        siteModelList = siteModelListInput.get(this);

                        alphaList.addParameter(alphaProposal);
                        invPrList.addParameter(invPrProposal);
                        ratesList.addParameter(ratesProposal);
                        siteModelList.addParameter(siteModelProposal);

                        //System.out.println("add "+ (++tempCount));
                        //System.out.println("sampler: "+ratesList.getID()+" "+ ratesList.getDimension());
                       //System.out.println("add "+ ratesList.getDimension());
                    }

                }

                //If any cluster has no member then it is removed.
                if(zeroCount > -1){

                    alphaList = alphaListInput.get(this);
                    invPrList = invPrListInput.get(this);
                    ratesList = ratesListInput.get(this);
                    siteModelList = siteModelListInput.get(this);

                    alphaList.removeParameter(zeroCount);
                    invPrList.removeParameter(zeroCount);
                    ratesList.removeParameter(zeroCount);
                    siteModelList.removeParameter(zeroCount);
                    //System.out.println("remove "+ (--tempCount));
                    //System.out.println("sampler: "+ratesList.getID()+" "+ ratesList.getDimension());

                }

            } /*else{
                System.out.println("No change");
            }*/




        }catch(Exception e){
                throw new RuntimeException(e);
            }
        return Double.POSITIVE_INFINITY;
    }

    public double getSiteLogLikelihood(
            double alpha,
            double invPr,
            double rates,
            double siteModel,
            int clusterID,
            int siteIndex,
            GeneralUnitSepTempTreeLikelihood tempLik){


            //System.out.println("flag 2");
            double siteLogLik =  dpTreeLikelihood.getSiteLogLikelihood(
                    DPNtdRateSepSiteModel.RATES,
                    clusterID,
                    siteIndex);
            if(Double.isNaN(siteLogLik)){
                //System.out.println("flag 2");
                siteLogLik =  tempLik.calculateLogP(
                        alpha,
                        invPr,
                        rates,
                        siteModel,
                        siteIndex
                );//todo need a new temp likelihood

            }

        return siteLogLik;
    }

    public QuietRealParameter[] getSamples(
            ParametricDistribution distr,
            double upper,
            double lower) throws Exception{
        QuietRealParameter[] samples = new QuietRealParameter[sampleSize];
        Double[][] sampleVals = distr.sample(sampleSize);
        for(int i = 0; i < samples.length;i++){
            samples[i] = new QuietRealParameter(sampleVals[i]);
            samples[i].setUpper(upper);
            samples[i].setLower(lower);
        }
        return samples;
    }

    @Override
	public void init(PrintStream out) throws Exception {
        out.print("count(Add cluster)\t");
    }

    @Override
	public void log(int nSample, PrintStream out) {
    	out.print((double)addClusterCount/(double)nSample+ "\t");
	}


    @Override
	public void close(PrintStream out) {
	}


}
