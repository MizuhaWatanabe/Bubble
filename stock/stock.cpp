#include <iostream>
#include <list>
#include <locale.h>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <curl/curl.h>

#define URL_BASE "http://stocks.finance.yahoo.co.jp/stocks/qi/?ids="
#define MATCH_1 "<td class=\"center yjM\">"
#define MATCH_2 "<strong class=\"yjMt\">"
#define MATCH_3 "<span class=\"yjSt profile\">"
#define REG_1 ">([^<>\n]+)<"

using namespace std;

/**
 * split関数
 * @param string str 分割したい文字列
 * @param string delim デリミタ
 * @return list<string> 分割された文字列
 */
list<string> split(string str, string delim) {
    list<string> result;
    int cutAt;
    while( (cutAt = str.find_first_of(delim)) != str.npos )
    {
        if(cutAt > 0)
        {
            result.push_back(str.substr(0, cutAt));
        }
        str = str.substr(cutAt + 1);
    }
    if(str.length() > 0)
    {
        result.push_back(str);
    }
    return result;
}

/*
 * CURL コールバック関数
 */
size_t callBackFunc(char* ptr, size_t size, size_t nmemb, string* stream) {
    int realsize = size * nmemb;
    stream->append(ptr, realsize);
    return realsize;
}

/*
 * [ CLASS ] Web Data
 */
class Web {

    // 各種宣言
    int intPage;                  // ページ番号
    char strPage[4];              // ページ番号(URL用)
    int intCountPage;             // ページ内件数
    int intCountCat;              // 業種内件数
    int intCountAll;              // 総件数
    bool blnFlagEnd;              // LOOP終了フラグ
    const char* pIter;            // イテレータ用ポインタ
    const char* pStockCode;       // 銘柄コードポインタ
    const char* pMarketName;      // 市場名ポインタ
    const char* pStockName;       // 銘柄名ポインタ
    const char* pCharacteristic;  // 特色ポインタ
    char strStockCode[5];         // 銘柄コード
    char strMarketName[32];       // 市場名
    char strStockName[128];       // 銘柄名
    char strCharacteristic[512];  // 特色
    int intI;                     // LOOP処理用インデックス

    // 業種コード定義
    const char *aryCategory[33] = {
        "0050", "1050", "2050", "3050", "3100",
        "3150", "3200", "3250", "3300", "3350",
        "3400", "3450", "3500", "3550", "3600",
        "3650", "3700", "3750", "3800", "4050",
        "5050", "5100", "5150", "5200", "5250",
        "6050", "6100", "7050", "7100", "7150",
        "7200", "8050", "9050"
    };

    public:

        // Constructor
        Web();

        // Webデータ読込
        int readHtml();

};

/*
 * [ Web ] Constructor
 */
Web::Web() {
    intCountAll = 0;
}

/*
 * [ Web ] HTML読込
 */
int Web::readHtml() {

    // 全業種コード分LOOP
    for (int i = 0; i < sizeof(aryCategory) / sizeof(aryCategory[0]); i++) {

        // ページ番号初期化
        intPage = 1;
        // 業種内件数初期化
        intCountCat = 0;
        // 終了フラグ初期化
        blnFlagEnd = false;

        while (blnFlagEnd == false) {

            // ページ内件数初期化
            intCountPage = 0;

            // URL 組み立て
            char strUrl[100] = "";
            sprintf(strPage, "%d", intPage);
            strcat(strUrl, URL_BASE);
            strcat(strUrl, aryCategory[i]);
            strcat(strUrl, "&p=");
            strcat(strUrl, strPage);

            // CURL 用宣言
            CURL *curl;
            CURLcode res;
            string chunk;

            // CURL 初期処理
            curl = curl_easy_init();
            // CURLOPT_URL に URL を指定
            curl_easy_setopt(curl, CURLOPT_URL, strUrl);
            // CURLOPT_WRITEFUNCTION にコールバック関数を指定
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, callBackFunc);
            // CURLOPT_WRITEDATA にコールバック関数にて処理されたあとのデータ格納ポインタを指定
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, (string*)&chunk);
            // リクエスト実行
            res = curl_easy_perform(curl);
            // リクエスト処理 OK 時
            if (res == CURLE_OK) {
                // split の実行
                list<string> strList = split(chunk, "\n");
                // イテレータの取得
                list<string>::iterator iter = strList.begin();
                // コンソール出力
                while( iter != strList.end() ) // 最後まで
                {
                    string strIter = *iter;
                    if (strIter.find(MATCH_1) == 0) {

                        // ページインクリメント
                        intCountPage++;
                        // イテレータのポインタ
                        pIter = strIter.c_str();

                        // 銘柄コード
                        pStockCode = pIter + 59;
                        for (int intI = 0; intI < 4; intI++) {
                            strStockCode[intI] = pStockCode[intI];
                        }
                        strStockCode[4] = '\0';

                        // 市場名
                        pMarketName = pStockCode + 37;
                        intI = 0;
                        while (intI != -1) {
                            if (pMarketName[intI] == '<') {
                                strMarketName[intI] = '\0';
                                intI = -1;
                            } else {
                                strMarketName[intI] = pMarketName[intI];
                                intI++;
                            }
                        }

                        // 銘柄名
                        pStockName = strstr(pMarketName, MATCH_2) + 57;
                        intI = 0;
                        while (intI != -1) {
                            if (pStockName[intI] == '<') {
                                strStockName[intI] = '\0';
                                intI = -1;
                            } else {
                                strStockName[intI] = pStockName[intI];
                                intI++;
                            }
                        }
                    }

                    if (strIter.find(MATCH_3) == 0) {

                        // イテレータのポインタ
                        pIter = strIter.c_str();

                        // 特色
                        pCharacteristic = pIter + 27;
                        intI = 0;
                        while (intI != -1) {
                            if (pCharacteristic[intI] == '<') {
                                strCharacteristic[intI] = '\0';
                                intI = -1;
                            } else {
                                strCharacteristic[intI] = pCharacteristic[intI];
                                intI++;
                            }
                        }
                        cout << "[" << aryCategory[i] << "] " << strStockCode << "[" << strMarketName << "] "
                             << strStockName << " : " << strCharacteristic << endl;
                    }

                    ++iter;
                }
                // 終了判定
                if (intCountPage == 20) {
                    intPage++;
                } else {
                    blnFlagEnd = true;
                }

            // リクエスト処理 NG 時
            } else {
                cout << "CURL ERROR!" << endl;
                return 1;
            }

            // ハンドラのクリーンアップ
            curl_easy_cleanup(curl);

            // LOOP 終了判定
            if (intCountPage == 0) {
                // 終了フラグ
                blnFlagEnd = true;
            } else {
                // 業種内件数
                intCountCat += intCountPage;
            }

        }

        // 総取得件数
        intCountAll += intCountCat;

    }

    // 総取得件数リターン
    return intCountAll;
}

/*
 * メイン処理
 */
int main()
{

    try {

        cout << "\n====< START >====\n" << endl;

        // Webデータ読込
        Web objWeb;
        int intCountAll = objWeb.readHtml();

        // 結果出力
        cout << "読込件数 : " << intCountAll << endl;

        cout << "\n====< E N D >====\n" << endl;

    }
    catch ( runtime_error ex ) {
        cout << ex.what() << endl;
    }

    return 0;

}