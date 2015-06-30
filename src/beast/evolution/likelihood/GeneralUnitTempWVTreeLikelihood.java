package beast.evolution.likelihood;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.GeneralUnitAlignment;
import beast.evolution.sitemodel.DummySiteModel;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 15/07/13
 * Time: 7:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneralUnitTempWVTreeLikelihood extends ExtendedTempWVTreeLikelihood {

    public void initAndValidate() throws Exception{

        super.initAndValidate();
        if(!(data instanceof GeneralUnitAlignment)){
            throw new RuntimeException("Data has to be an instance of GeneralUnitAlignment object.");
        }
    }


    public double[] calculateLogP (
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            RealParameter rate,
            int[] units) throws Exception{
        GeneralUnitAlignment alignment  = (GeneralUnitAlignment)data;


        double[] siteLogP = new double[units.length];
        calculateLogP(modelParameters,modelCode,freqs,rate);
        for(int i = 0; i < units.length; i++){
            int[] sites = alignment.getSitesByUnit(units[i]);
            for(int site:sites){
                siteLogP[i] += patternLogLikelihoods[dataInput.get().getPatternIndex(site)];
            }
        }
        return siteLogP;
    }

    public double[] calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            RealParameter rate,
            int[] units,
            int exceptUnit) {
            double[] siteLogP = new double[units.length-1];
        GeneralUnitAlignment alignment  = (GeneralUnitAlignment)data;
        try{

            calculateLogP(modelParameters,modelCode,freqs,rate);
            int k = 0;
            for(int i = 0; i < units.length;i++){
                if(units[i] != exceptUnit){
                    //System.out.println(sites[i] +" "+exceptSite);
                    int[] sites = alignment.getSitesByUnit(units[i]);
                    for(int site:sites){
                        siteLogP[k] += patternLogLikelihoods[dataInput.get().getPatternIndex(site)];
                    }
                    k++;
                }

            }


        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return siteLogP;
    }

    public void setupPatternWeightsFromSites(int[] units){
        int[] tempWeights = new int[dataInput.get().getPatternCount()];
        int patIndex;

        for(int i = 0; i < units.length; i++){
            int[] sites = ((GeneralUnitAlignment)data).getSitesByUnit(units[i]);
            for(int j = 0; j < sites.length; j++){
                patIndex = dataInput.get().getPatternIndex(sites[j]);
                tempWeights[patIndex] = 1;
            }
        }

        setPatternWeights(tempWeights);
    }



    public double[] calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int[] units,
            int exceptUnit) throws Exception{
            double[] siteLogP = new double[units.length-1];
        try{

            GeneralUnitAlignment alignment = (GeneralUnitAlignment)data;

            calculateLogP(modelParameters,modelCode,freqs);
            int k = 0;
            for(int i = 0; i < units.length;i++){
                if(units[i] != exceptUnit){
                    int[] sites = alignment.getSitesByUnit(units[i]);
                    for(int site:sites){
                        siteLogP[k] += patternLogLikelihoods[dataInput.get().getPatternIndex(site)];
                    }
                    k++;
                }
            }


        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return siteLogP;
    }


    public double[] calculateLogP(
            RealParameter rateParameter,
            int[] units) throws Exception{
        double[] siteLogP = new double[units.length];
        GeneralUnitAlignment alignment = (GeneralUnitAlignment)data;
        int k = 0;
        try{
            ((DummySiteModel)m_siteModel).getRateParameter().setValueQuietly(0,rateParameter.getValue());
            calculateLogP();
            for(int i = 0; i < units.length;i++){
                int[] sites = alignment.getSitesByUnit(units[i]);
                for(int site:sites){
                    siteLogP[k] += patternLogLikelihoods[dataInput.get().getPatternIndex(site)];
                }
                k++;
                //System.out.println(siteLogP[i]);
            }
        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return siteLogP;
    }

    public double[] calculateLogP(
            RealParameter rateParameter,
            int[] units,
            int exceptUnit) throws Exception{
        GeneralUnitAlignment alignment = (GeneralUnitAlignment)data;
        double[] siteLogP = new double[units.length-1];
        try{
            ((DummySiteModel)m_siteModel).getRateParameter().setValueQuietly(0,rateParameter.getValue());
            calculateLogP();
            int k = 0;
            for(int i = 0; i < units.length;i++){
                if(units[i] != exceptUnit){
                    int[] sites = alignment.getSitesByUnit(units[i]);
                    for(int site:sites){
                        siteLogP[k] += patternLogLikelihoods[dataInput.get().getPatternIndex(site)];
                    }
                    k++;
                }
            }
        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return siteLogP;
    }


    public double[] calculateLogP (
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            RealParameter alpha,
            RealParameter invPr,
            RealParameter rate,
            RealParameter siteModelChoice,
            int[] units) throws Exception{
        GeneralUnitAlignment alignment  = (GeneralUnitAlignment)data;


        double[] siteLogP = new double[units.length];
        calculateLogP(
                modelParameters,
                modelCode,
                freqs,
                alpha,
                invPr,
                rate,
                siteModelChoice
        );
        for(int i = 0; i < units.length; i++){
            int[] sites = alignment.getSitesByUnit(units[i]);
            for(int site:sites){
                siteLogP[i] += patternLogLikelihoods[dataInput.get().getPatternIndex(site)];
            }
        }
        return siteLogP;
    }

    public double[] calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            RealParameter alpha,
            RealParameter invPr,
            RealParameter rate,
            RealParameter siteModelChoice,
            int[] units,
            int exceptUnit) {
            double[] siteLogP = new double[units.length-1];
        GeneralUnitAlignment alignment  = (GeneralUnitAlignment)data;
        try{

            calculateLogP(
                    modelParameters,
                    modelCode,
                    freqs,
                    alpha,
                    invPr,
                    rate,
                    siteModelChoice
            );

            int k = 0;
            for(int i = 0; i < units.length;i++){
                if(units[i] != exceptUnit){
                    //System.out.println(sites[i] +" "+exceptSite);
                    int[] sites = alignment.getSitesByUnit(units[i]);
                    for(int site:sites){
                        siteLogP[k] += patternLogLikelihoods[dataInput.get().getPatternIndex(site)];
                    }
                    k++;
                }

            }


        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return siteLogP;
    }






}
