package edu.albany.cs.base;

/**
 * Created by baojian bzhou6@albany.edu on 2/4/17.
 */

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentMap;

import org.python.core.PyDictionary;
import org.python.core.PyFile;
import org.python.core.PyList;
import org.python.core.PyObject;
import org.python.core.PyString;
import org.python.modules.cPickle;

public class PickleLoader {

    private final String filePath;

    public PickleLoader() {
        filePath = null;
    }

    public PickleLoader(String filePath) {
        this.filePath = filePath;
        loadData();
    }

    public PickleLoader(Path filePath) {
        this(filePath.toString());
    }

    private void loadData() {
        HashMap<String, ArrayList<String>> idToCountries = new HashMap<>();
        File f = new File(filePath);
        BufferedReader bufR;
        StringBuilder strBuilder = new StringBuilder();
        try {
            bufR = new BufferedReader(new FileReader(f));
            String aLine;
            while (null != (aLine = bufR.readLine())) {
                strBuilder.append(aLine).append("n");
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        PyString pyStr = new PyString(strBuilder.toString());
        PyDictionary idToCountriesObj = (PyDictionary) cPickle.loads(pyStr);
        ConcurrentMap<PyObject, PyObject> aMap = idToCountriesObj.getMap();
        for (Map.Entry<PyObject, PyObject> entry : aMap.entrySet()) {
            String appId = entry.getKey().toString();
            PyList countryIdList = (PyList) entry.getValue();
            List<String> countryList = (List<String>) countryIdList.subList(0, countryIdList.size());
            ArrayList<String> countryArrList = new ArrayList<String>(countryList);
            idToCountries.put(appId, countryArrList);
        }
        Utils.stop();
    }

    public static void main(String[] args) {
        new PickleLoader("/home/baojian/Dropbox/expriments/ACES_Data/data/U133A_combat_DMFS_nwEdgesKEGG_data.pkl");
    }


    public HashMap<String, ArrayList<String>> getIdToCountriesFileStream() {
        HashMap<String, ArrayList<String>> idToCountries = new HashMap<String, ArrayList<String>>();
        File f = new File("test.pkl");
        InputStream fs = null;
        try {
            fs = new FileInputStream(f);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return null;
        }
        PyFile pyStr = new PyFile(fs);
        PyDictionary idToCountriesObj = (PyDictionary) cPickle.load(pyStr);
        ConcurrentMap<PyObject, PyObject> aMap = idToCountriesObj.getMap();
        for (Map.Entry<PyObject, PyObject> entry : aMap.entrySet()) {
            String appId = entry.getKey().toString();
            PyList countryIdList = (PyList) entry.getValue();
            List<String> countryList = (List<String>) countryIdList.subList(0, countryIdList.size());
            ArrayList<String> countryArrList = new ArrayList<String>(countryList);
            idToCountries.put(appId, countryArrList);
        }
//        System.out.println(idToCountries.toString());
        return idToCountries;
    }
}