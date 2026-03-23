//! AI integration — daimon/hoosh client for prakash.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DaimonConfig {
    pub endpoint: String,
    pub api_key: Option<String>,
}

impl Default for DaimonConfig {
    fn default() -> Self {
        Self { endpoint: "http://localhost:8090".into(), api_key: None }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HooshConfig {
    pub endpoint: String,
}

impl Default for HooshConfig {
    fn default() -> Self {
        Self { endpoint: "http://localhost:8088".into() }
    }
}

pub struct DaimonClient {
    config: DaimonConfig,
    client: reqwest::Client,
}

impl DaimonClient {
    pub fn new(config: DaimonConfig) -> Self {
        Self {
            config,
            client: reqwest::Client::builder()
                .timeout(std::time::Duration::from_secs(30))
                .build()
                .expect("failed to build HTTP client"),
        }
    }

    pub async fn register_agent(&self) -> anyhow::Result<String> {
        let body = serde_json::json!({
            "name": "prakash",
            "capabilities": ["optics", "ray_tracing", "spectral", "pbr", "wave_optics"],
        });
        let resp = self.client
            .post(format!("{}/v1/agents/register", self.config.endpoint))
            .json(&body)
            .send()
            .await?;
        let data: serde_json::Value = resp.json().await?;
        Ok(data["agent_id"].as_str().unwrap_or("unknown").to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let c = DaimonConfig::default();
        assert_eq!(c.endpoint, "http://localhost:8090");
    }

    #[test]
    fn test_hoosh_default() {
        let c = HooshConfig::default();
        assert_eq!(c.endpoint, "http://localhost:8088");
    }
}
