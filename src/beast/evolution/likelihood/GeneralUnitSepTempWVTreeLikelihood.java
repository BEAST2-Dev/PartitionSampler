package beast.evolution.likelihood;

import beast.core.parameter.QuietRealParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.GeneralUnitAlignment;
import beast.evolution.sitemodel.DPNtdBMAGIBMASepSiteModel;
import beast.evolution.sitemodel.DPNtdRateSepSiteModel;
import beast.evolution.sitemodel.DummySiteModel;
import beast.evolution.sitemodel.QuietGammaSiteBMA;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 17/07/13
 * Time: 3:36 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneralUnitSepTempWVTreeLikelihood extends ExtendedSepTempWVTreeLikelihood {

    public double[] calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int[] unitIndexP,
            int exceptUnit){

        int[] units = new int[unitIndexP.length - 1];
        int k = 0;
        for(int i = 0; i < unitIndexP.length;i++){
            if(unitIndexP[i] != exceptUnit){
                units[k++] = unitIndexP[i];
            }
        }
        return calculateLogP(modelParameters,modelCode,freqs,units);

    }

    public double[] calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int[] units){

        if(dpNtdRateSepSiteModel instanceof DPNtdBMAGIBMASepSiteModel){

            return calculateLogP2(modelParameters,modelCode,freqs,units);

        }else{
            return calculateLogP1(modelParameters,modelCode,freqs,units);
        }

    }

    public double[] calculateLogP1(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int[] units){

        HashMap<Integer, Integer> unitOrderMap = new HashMap<Integer, Integer>();
        for(int i = 0; i < units.length; i++){
            unitOrderMap.put(units[i],i);
        }

        GeneralUnitAlignment alignment = (GeneralUnitAlignment) dataInput.get();
        double[] logPs = new double[units.length];
        setModelParameterVals(modelParameters,modelCode,freqs);

        try{
            //Stores a list of rates
            ArrayList<QuietRealParameter> rates = new ArrayList<QuietRealParameter>();
            //Stores the ID number of the substModel and pattern weights
            HashMap<Integer,int[]> clusterWeightsMap = new HashMap<Integer,int[]>();
            //Stores the ID number of the substModel and the sites with that model
            HashMap<Integer, ArrayList<Integer>> clusterSitesMap = new HashMap<Integer, ArrayList<Integer>>();

            int patternIndex;
            int rateID;

            for(int i = 0; i < units.length;i++){
                int[] sites = alignment.getSitesByUnit(units[i]);



                QuietRealParameter rate = dpNtdRateSepSiteModel.getRate(units[i]);
                rateID = rate.getIDNumber();


                for(int site:sites){


                    patternIndex = alignment.getPatternIndex(site);

                    if(clusterWeightsMap.containsKey(rateID)){
                        //increment the weight of the pattern at site siteIndex[i] by 1
                        clusterWeightsMap.get(rateID)[patternIndex]++;
                        //add the site index to the cluster
                        clusterSitesMap.get(rateID).add(site);

                    }else{
                        //Create a new pattern weights vector
                        int[] patternWeights = new int[alignment.getPatternCount()];
                        //Increment the weight of the pattern at site siteIndex[i]
                        patternWeights[patternIndex]++;
                        //Add the new pattern weight vector to the HashMap
                        clusterWeightsMap.put(rateID,patternWeights);

                        //Add the new substModel to the array list
                        rates.add(rate);

                        //Create a new array list to store sites
                        ArrayList<Integer> clusterSites = new ArrayList<Integer>();
                        //Add site index to the corresponding list
                        clusterSites.add(site);
                        //Put the list in the HashMap
                        clusterSitesMap.put(rateID,clusterSites);
                    }
                }

            }

            int rateCount = rates.size();
            int unit;

            for(int i = 0; i < rateCount;i++){

                m_siteModel.substModelInput.setValue(substModel, m_siteModel);
                ((DummySiteModel)m_siteModel).getRateParameter().setValueQuietly(0,rates.get(i).getValue());
                rateID = rates.get(i).getIDNumber();
                setPatternWeights(clusterWeightsMap.get(rateID));
                calculateLogP();

                //Get all the site indices in this rate category
                ArrayList<Integer> sites = clusterSitesMap.get(rateID);
                for(Integer site:sites){
                    unit = alignment.getUnitBySite(site);
                    logPs[unitOrderMap.get(unit)] += patternLogLikelihoods[dataInput.get().getPatternIndex(site)] ;
                }
            }

        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return logPs;
    }


    public double[] calculateLogP2(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int[] units){



        HashMap<Integer, Integer> unitOrderMap = new HashMap<Integer, Integer>();
        for(int i = 0; i < units.length; i++){
            //System.out.println("unit: "+units[i]);
            unitOrderMap.put(units[i],i);
        }

        GeneralUnitAlignment alignment = (GeneralUnitAlignment) dataInput.get();
        double[] logPs = new double[units.length];
        setModelParameterVals(modelParameters,modelCode,freqs);

        try{
            //Stores a list of rates
            ArrayList<Double> rates = new ArrayList<Double>();
            ArrayList<Double> alpha = new ArrayList<Double>();
            ArrayList<Double> invPr = new ArrayList<Double>();
            ArrayList<Double> siteModelChoice = new ArrayList<Double>();
            ArrayList<Integer> rateIDs = new ArrayList<Integer>();

            DPNtdBMAGIBMASepSiteModel dpSiteModel = (DPNtdBMAGIBMASepSiteModel) dpNtdRateSepSiteModel;


            //Stores the ID number of the substModel and pattern weights
            HashMap<Integer,int[]> clusterWeightsMap = new HashMap<Integer,int[]>();
            //Stores the ID number of the substModel and the sites with that model
            HashMap<Integer, ArrayList<Integer>> clusterSitesMap = new HashMap<Integer, ArrayList<Integer>>();

            int patternIndex;
            int rateID;

            for(int i = 0; i < units.length;i++){
                //System.out.println("flag1, unit: "+units[i]);
                int[] sites = alignment.getSitesByUnit(units[i]);



                //QuietRealParameter rate = dpNtdRateSepSiteModel.getRate(units[i]);
                double rateVal = dpSiteModel.getRateValue(units[i]);
                double alphaVal = dpSiteModel.getAlphaValue(units[i]);
                double invPrVal = dpSiteModel.getInvPrValue(units[i]);
                double siteModelChoiceVal = dpSiteModel.getModelChoiceValue(units[i]);
                //rateID = rate.getIDNumber();
                rateID = dpSiteModel.getRateID(units[i]);
                //System.out.println("rate ID: "+rateID);

                for(int site:sites){
                    //System.out.println("flag1, site: "+site);


                    patternIndex = alignment.getPatternIndex(site);

                    if(clusterWeightsMap.containsKey(rateID)){
                        //increment the weight of the pattern at site siteIndex[i] by 1
                        clusterWeightsMap.get(rateID)[patternIndex]++;
                        //add the site index to the cluster
                        clusterSitesMap.get(rateID).add(site);

                    }else{
                        //Create a new pattern weights vector
                        int[] patternWeights = new int[alignment.getPatternCount()];
                        //Increment the weight of the pattern at site siteIndex[i]
                        patternWeights[patternIndex]++;
                        //Add the new pattern weight vector to the HashMap
                        clusterWeightsMap.put(rateID,patternWeights);

                        //Add the new substModel to the array list
                        rates.add(rateVal);
                        alpha.add(alphaVal);
                        invPr.add(invPrVal);
                        siteModelChoice.add(siteModelChoiceVal);
                        rateIDs.add(rateID);

                        //Create a new array list to store sites
                        ArrayList<Integer> clusterSites = new ArrayList<Integer>();
                        //Add site index to the corresponding list
                        clusterSites.add(site);
                        //Put the list in the HashMap
                        clusterSitesMap.put(rateID,clusterSites);
                    }
                }

            }

            int rateCount = rates.size();
            int unit;

            for(int i = 0; i < rateCount;i++){

                m_siteModel.substModelInput.setValue(substModel, m_siteModel);
                //((DummySiteModel)m_siteModel).getRateParameter().setValueQuietly(0,rates.get(i).getValue());
                ((QuietGammaSiteBMA)m_siteModel).setMuValueQuietly(rates.get(i));
                ((QuietGammaSiteBMA)m_siteModel).setShapeValueQuietly(alpha.get(i));
                ((QuietGammaSiteBMA)m_siteModel).setInvPrValueQuietly(invPr.get(i));
                ((QuietGammaSiteBMA)m_siteModel).setModelChoiceQuietly(siteModelChoice.get(i));
                ((QuietGammaSiteBMA)m_siteModel).setRatesKnown(false);




                //rateID = rates.get(i).getIDNumber();
                setPatternWeights(clusterWeightsMap.get(rateIDs.get(i)));
                calculateLogP();

                //Get all the site indices in this rate category
                ArrayList<Integer> sites = clusterSitesMap.get(rateIDs.get(i));
                for(Integer site:sites){
                    //System.out.println("site: "+site);
                    unit = alignment.getUnitBySite(site);
                    //System.out.println("unit: "+unit);
                    logPs[unitOrderMap.get(unit)] += patternLogLikelihoods[dataInput.get().getPatternIndex(site)] ;
                }
            }

        }catch(Exception e){
            throw new RuntimeException(e);

        }
        /*for(double logP:logPs){
            System.out.println("logP: "+logP);
        }*/
        return logPs;
    }

    public double[] calculateLogP(
            RealParameter rateParameter,
            int[] unitIndexP,
            int unitExcept)throws Exception{

        int[] units = new int[unitIndexP.length - 1];
        int k = 0;
        for(int i = 0; i < unitIndexP.length;i++){
            if(unitIndexP[i] != unitExcept){
                units[k++] =  unitIndexP[i];
            }
        }
        return calculateLogP(rateParameter,units);

    }


    public double[] calculateLogP (
            RealParameter rateParameter,
            int[] units) throws Exception{

        HashMap<Integer, Integer> unitOrderMap = new HashMap<Integer,Integer>();
        for(int i = 0;i < units.length;i++){
            //System.out.println("units: "+units[i]);
            unitOrderMap.put(units[i], i);
        }
        GeneralUnitAlignment alignment = (GeneralUnitAlignment) dataInput.get();

        double[] logPs = new double[units.length];



        //Stores a list of subst models
        ArrayList<SwitchingNtdBMA> ntdBMAMap = new ArrayList<SwitchingNtdBMA>();
        //Stores the ID number of the substModel and pattern weights
        HashMap<Integer,int[]> clusterWeightsMap = new HashMap<Integer,int[]>();
        //Stores the ID number of the substModel and the sites with that model
        HashMap<Integer, ArrayList<Integer>> clusterSitesMap = new HashMap<Integer, ArrayList<Integer>>();

        int patternIndex;
        int modelID;

        for(int i = 0; i < units.length;i++){

            //System.out.println((SwitchingNtdBMA)dpNtdRateSepSiteModel.getModel(siteIndex[i]));
            SwitchingNtdBMA ntdBMA = ((SwitchingNtdBMA)dpNtdRateSepSiteModel.getModel(units[i]));
            modelID = ntdBMA.getIDNumber();
            int[] sites = alignment.getSitesByUnit(units[i]);

            for(int site:sites){

                patternIndex = alignment.getPatternIndex(site);

                if(clusterWeightsMap.containsKey(modelID)){
                    //increment the weight of the pattern at jth site of unit[i] by 1
                    clusterWeightsMap.get(modelID)[patternIndex]++;
                    //add the site index in unit i to the cluster
                    clusterSitesMap.get(modelID).add(site);

                }else{
                    //Create a new pattern weights vector
                    int[] patternWeights = new int[alignment.getPatternCount()];
                    //Increment the weight of the pattern at jth site of unit[i] by 1
                    patternWeights[patternIndex]++;
                    //Add the new pattern weight vector to the HashMap
                    clusterWeightsMap.put(modelID,patternWeights);

                    //Add the new substModel to the array list
                    ntdBMAMap.add(ntdBMA);

                    //Create a new array list to store sites
                    ArrayList<Integer> clusterSites = new ArrayList<Integer>();
                    //Add site index in unit i to the corresponding list
                    clusterSites.add(site);
                    //Put the list in the HashMap
                    clusterSitesMap.put(modelID,clusterSites);
                }
            }

        }

        int substModelCount = ntdBMAMap.size();
        int unit;

        for(int i = 0; i < substModelCount;i++){

            m_siteModel.substModelInput.setValue(ntdBMAMap.get(i), m_siteModel);
            ((DummySiteModel)m_siteModel).getRateParameter().setValueQuietly(0,rateParameter.getValue());
            modelID = ntdBMAMap.get(i).getIDNumber();
            setPatternWeights(clusterWeightsMap.get(modelID));
            calculateLogP();

            ArrayList<Integer> sites = clusterSitesMap.get(modelID);
            for(Integer site:sites){
                unit = alignment.getUnitBySite(site);
                //System.out.println("unit: "+unit);
                logPs[unitOrderMap.get(unit)] += patternLogLikelihoods[dataInput.get().getPatternIndex(site)] ;
            }
        }
        //System.out.println(logP);


        return logPs;

    }


    public double[] calculateLogP(
            double alpha,
            double invPr,
            double rate,
            double siteModelChoice,
            int[] unitIndexP,
            int unitExcept) throws Exception{

        int[] units = new int[unitIndexP.length - 1];
        int k = 0;
        for(int i = 0; i < unitIndexP.length;i++){
            if(unitIndexP[i] != unitExcept){
                units[k++] =  unitIndexP[i];
            }
        }
        return calculateLogP(
                alpha,
                invPr,
                rate,
                siteModelChoice,
                units
        );


    }




    public double[] calculateLogP (
            double alpha,
            double invPr,
            double rate,
            double siteModelChoice,
            int[] units) throws Exception{

        HashMap<Integer, Integer> unitOrderMap = new HashMap<Integer,Integer>();
        for(int i = 0;i < units.length;i++){
            //System.out.println("units: "+units[i]);
            unitOrderMap.put(units[i], i);
        }
        GeneralUnitAlignment alignment = (GeneralUnitAlignment) dataInput.get();

        double[] logPs = new double[units.length];



        //Stores a list of subst models
        ArrayList<SwitchingNtdBMA> ntdBMAMap = new ArrayList<SwitchingNtdBMA>();
        //Stores the ID number of the substModel and pattern weights
        HashMap<Integer,int[]> clusterWeightsMap = new HashMap<Integer,int[]>();
        //Stores the ID number of the substModel and the sites with that model
        HashMap<Integer, ArrayList<Integer>> clusterSitesMap = new HashMap<Integer, ArrayList<Integer>>();

        int patternIndex;
        int modelID;

        for(int i = 0; i < units.length;i++){

            //System.out.println((SwitchingNtdBMA)dpNtdRateSepSiteModel.getModel(siteIndex[i]));
            SwitchingNtdBMA ntdBMA = ((SwitchingNtdBMA)dpNtdRateSepSiteModel.getModel(units[i]));
            modelID = ntdBMA.getIDNumber();
            int[] sites = alignment.getSitesByUnit(units[i]);

            for(int site:sites){

                patternIndex = alignment.getPatternIndex(site);

                if(clusterWeightsMap.containsKey(modelID)){
                    //increment the weight of the pattern at jth site of unit[i] by 1
                    clusterWeightsMap.get(modelID)[patternIndex]++;
                    //add the site index in unit i to the cluster
                    clusterSitesMap.get(modelID).add(site);

                }else{
                    //Create a new pattern weights vector
                    int[] patternWeights = new int[alignment.getPatternCount()];
                    //Increment the weight of the pattern at jth site of unit[i] by 1
                    patternWeights[patternIndex]++;
                    //Add the new pattern weight vector to the HashMap
                    clusterWeightsMap.put(modelID,patternWeights);

                    //Add the new substModel to the array list
                    ntdBMAMap.add(ntdBMA);

                    //Create a new array list to store sites
                    ArrayList<Integer> clusterSites = new ArrayList<Integer>();
                    //Add site index in unit i to the corresponding list
                    clusterSites.add(site);
                    //Put the list in the HashMap
                    clusterSitesMap.put(modelID,clusterSites);
                }
            }

        }

        int substModelCount = ntdBMAMap.size();
        int unit;

        for(int i = 0; i < substModelCount;i++){

            m_siteModel.substModelInput.setValue(ntdBMAMap.get(i), m_siteModel);
            //((DummySiteModel)m_siteModel).getRateParameter().setValueQuietly(0,rateParameter.getValue());
            ((QuietGammaSiteBMA)m_siteModel).setMuValueQuietly(rate);
            ((QuietGammaSiteBMA)m_siteModel).setShapeValueQuietly(alpha);
            ((QuietGammaSiteBMA)m_siteModel).setInvPrValueQuietly(invPr);
            ((QuietGammaSiteBMA)m_siteModel).setModelChoiceQuietly(siteModelChoice);
            ((QuietGammaSiteBMA)m_siteModel).setRatesKnown(false);

            modelID = ntdBMAMap.get(i).getIDNumber();
            setPatternWeights(clusterWeightsMap.get(modelID));
            calculateLogP();

            ArrayList<Integer> sites = clusterSitesMap.get(modelID);
            for(Integer site:sites){
                unit = alignment.getUnitBySite(site);
                //System.out.println("unit: "+unit);
                logPs[unitOrderMap.get(unit)] += patternLogLikelihoods[dataInput.get().getPatternIndex(site)] ;
            }
        }
        //System.out.println(logP);


        return logPs;

    }


}
